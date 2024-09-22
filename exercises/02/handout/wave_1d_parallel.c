#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>

// TASK: T1a
// Include the MPI headerfile
// BEGIN: T1a
#include <mpi.h>
// END: T1a


// Option to change numerical precision.
typedef int64_t int_t;
typedef double real_t;


// TASK: T1b
// Declare variables each MPI process will need
// BEGIN: T1b
int world_size, world_rank;
// END: T1b

int_t local_N, local_start;

// Simulation parameters: size, step count, and how often to save the state.
const int_t
N = 65536,
max_iteration = 100000,
snapshot_freq = 500;

// Wave equation parameters, time step is derived from the space step.
const real_t
c  = 1.0,
dx = 1.0;
real_t
dt;

// Buffers for three time steps, indexed with 2 ghost points for the boundary.
real_t
*buffers[3] = { NULL, NULL, NULL };

#define U_prv(i) buffers[0][(i)+1]
#define U(i)     buffers[1][(i)+1]
#define U_nxt(i) buffers[2][(i)+1]


// Convert 'struct timeval' into seconds in double prec. floating point
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)


// TASK: T8
// Save the present time step in a numbered file under 'data/'.
void domain_save ( int_t step )
{
    // BEGIN: T8
    if (world_rank == 0) {
        char filename[256];
        sprintf ( filename, "data/%.5ld.dat", step );
        FILE *out = fopen ( filename, "wb" );
        fwrite ( &U(0), sizeof(real_t), N, out );
        fclose ( out );
    }
    // END: T8
}


// TASK: T3
// Allocate space for each process' sub-grids
// Set up our three buffers, fill two with an initial cosine wave,
// and set the time step.
void domain_initialize ( void )
{
    // BEGIN: T3

    // printf("Rank %d: Entering domain_initialize\n", world_rank);

    // Calculate the local domain size for each process
    local_N = N / world_size;
    if (world_rank < N % world_size) {
        local_N++;
    }

    // Calculate the starting index for this process
    if (world_rank < N % world_size) {
        local_start = world_rank * (N / world_size) + world_rank;
    } else {
        local_start = world_rank * (N / world_size) + N % world_size;
    }

    // printf("Rank %d: local_N = %ld, local_start = %ld\n", world_rank, local_N, local_start);

    // Allocate memory for local buffers (including ghost cells)
    for (int i = 0; i < 3; i++) {
        if(world_rank == 0) { // Rank 0 has to have the same buffer size as the global domain
            buffers[i] = malloc((N + 2) * sizeof(real_t));
        } else {
            buffers[i] = malloc((local_N + 2) * sizeof(real_t));
        }
        // printf("Rank %d: Allocated buffer %d at address %p\n", world_rank, i, (void*)buffers[i]);
    }

    // Initialize the local domain with the cosine wave
    for (int_t i = 0; i < local_N; i++) {
        real_t x = (local_start + i) / (real_t) N;
        U_prv(i) = U(i) = cos(M_PI * x);
    }

    // Set the time step for 1D case.
    dt = dx / c;

    // printf("Rank %d: Exiting domain_initialize\n", world_rank);

    // END: T3
}


// Return the memory to the OS.
void domain_finalize ( void )
{
    // printf("Rank %d: Entering domain_finalize\n", world_rank);

    free(buffers[0]);
    free(buffers[1]);
    free(buffers[2]);

    // printf("Rank %d: Exiting domain_finalize\n", world_rank);
}


// Rotate the time step buffers.
void move_buffer_window ( void )
{
    real_t *temp = buffers[0];
    buffers[0] = buffers[1];
    buffers[1] = buffers[2];
    buffers[2] = temp; 
}


// TASK: T4
// Derive step t+1 from steps t and t-1.
void time_step ( void )
{
    // BEGIN: T4
    for ( int_t i = 0; i < local_N; i++ )
    {
        U_nxt(i) = -U_prv(i) + 2.0*U(i)
            + (dt*dt*c*c)/(dx*dx) * (U(i-1)+U(i+1)-2.0*U(i));
    }
    // END: T4
}


// TASK: T6
// Neumann (reflective) boundary condition.
void boundary_condition ( void )
{
    // BEGIN: T6
    if (world_rank == 0) {
        U(-1) = U(1);
    }
    if (world_rank == world_size - 1) {
        U(local_N) = U(local_N-2);
    }
    // END: T6
}



// TASK: T5
// Communicate the border between processes.
void border_exchange( void )
{
    // BEGIN: T5

    // Exchange data with the right neighbour
    if (world_rank < world_size - 1) {
        MPI_Send(&U(local_N-1), 1, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&U(local_N), 1, MPI_DOUBLE, world_rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Exchange data with the left neighbour
    if (world_rank > 0) {
        MPI_Send(&U(0), 1, MPI_DOUBLE, world_rank - 1, 1, MPI_COMM_WORLD);
        MPI_Recv(&U(-1), 1, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // END: T5
}



// TASK: T7
// Every process needs to communicate its results
// to root and assemble it in the root buffer
void send_data_to_root()
{
    // BEGIN: T7
    if (world_rank == 0) {

        real_t *complete_domain = malloc(N * sizeof(real_t));

        memcpy(complete_domain, &U(0), local_N * sizeof(real_t));

        int offset = local_N;

        for (int i = 1; i < world_size; i++) {
            int recv_size;

            // First, receive the size of data from the process
            MPI_Recv(&recv_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Then, receive the actual data
            MPI_Recv(complete_domain + offset, recv_size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            offset += recv_size;
        }

        // printf("Rank %d: Received data from all processes\n", world_rank);

        memcpy(&U(0), complete_domain, N * sizeof(real_t));

        // printf("Rank %d: Copied data to U buffer\n", world_rank);

        free(complete_domain);

        // printf("Rank %d: Freed global_U\n", world_rank);

    } else {
        // Send the size of the data to the root process
        MPI_Send(&local_N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        // Send the actual data to the root process
        MPI_Send(&U(0), local_N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    // END: T7

    // printf ( "Rank %d: Exiting send_data_to_root\n", world_rank );
}


// Main time integration.
void simulate( void )
{
    // Go through each time step.
    for ( int_t iteration=0; iteration<=max_iteration; iteration++ )
    {
        if ( (iteration % snapshot_freq)==0 )
        {
            // printf ( "Rank %d: Iteration %ld, calling send_data_to_root\n", world_rank, iteration );
            send_data_to_root();
            domain_save ( iteration / snapshot_freq );
        }

        // Derive step t+1 from steps t and t-1.
        border_exchange();
        boundary_condition();
        time_step();

        move_buffer_window();
    }
}


int main ( int argc, char **argv )
{
    // TASK: T1c
    // Initialise MPI
    // BEGIN: T1c
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // END: T1c

    struct timeval t_start, t_end;
    domain_initialize();

    // TASK: T2
    // Time your code
    // BEGIN: T2
    if (world_rank == 0) {
        gettimeofday(&t_start, NULL);
    }

    simulate();

    if (world_rank == 0) {
        gettimeofday(&t_end, NULL);
        double elapsed = WALLTIME(t_end) - WALLTIME(t_start);
        printf("Elapsed time: %f seconds\n", elapsed);
    }
    // END: T2

    domain_finalize();

    // TASK: T1d
    // Finalise MPI
    // BEGIN: T1d
    MPI_Finalize();
    // END: T1d

    exit ( EXIT_SUCCESS );
}

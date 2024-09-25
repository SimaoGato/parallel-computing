#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>

#include "argument_utils.h"

// TASK: T1a
// Include the MPI headerfile
#include <mpi.h>

// Convert 'struct timeval' into seconds in double prec. floating point
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

// Option to change numerical precision
typedef int64_t int_t;
typedef double real_t;

// Buffers for three time steps, indexed with 2 ghost points for the boundary
real_t *buffers[3] = { NULL, NULL, NULL };

// TASK: T1b
// Declare variables each MPI process will need
int world_size, world_rank;
int cart_dims[2], cart_periods[2] = {0, 0}, cart_coords[2];
MPI_Comm cart_comm;
int local_M, local_N;
int north, south, east, west;

#define U_prv(i,j) buffers[0][((i)+1)*(local_N+2)+(j)+1]
#define U(i,j)     buffers[1][((i)+1)*(local_N+2)+(j)+1]
#define U_nxt(i,j) buffers[2][((i)+1)*(local_N+2)+(j)+1]

// Simulation parameters: size, step count, and how often to save the state
int_t
    M = 256,    // rows
    N = 256,    // cols
    max_iteration = 4000,
    snapshot_freq = 20;

// Wave equation parameters, time step is derived from the space step
const real_t
    c  = 1.0,
    dx = 1.0,
    dy = 1.0;
real_t dt;

// Rotate the time step buffers.
void move_buffer_window(void)
{
    real_t *temp = buffers[0];
    buffers[0] = buffers[1];
    buffers[1] = buffers[2];
    buffers[2] = temp;
}

// TASK: T4
// Set up our three buffers, and fill two with an initial perturbation
// and set the time step.
void domain_initialize(void)
{
    local_M = M / cart_dims[0];
    local_N = N / cart_dims[1];

    // Initialize neighbour ranks
    north = south = east = west = MPI_PROC_NULL;
    MPI_Cart_shift(cart_comm, 0, 1, &north, &south);
    MPI_Cart_shift(cart_comm, 1, 1, &west, &east);

    // Print the domain decomposition
    if (world_rank == 0)
    {
        printf("Decomposition: %dx%d\n", cart_dims[0], cart_dims[1]);
    }

    printf("Rank %d: local_M = %d, local_N = %d\n", world_rank, local_M, local_N);
    printf("Rank %d: cart_coords = (%d, %d)\n", world_rank, cart_coords[0], cart_coords[1]);
    printf("Rank %d: north = %d, south = %d, east = %d, west = %d\n", world_rank, north, south, east, west);

    buffers[0] = malloc((local_M+2)*(local_N+2)*sizeof(real_t));
    buffers[1] = malloc((local_M+2)*(local_N+2)*sizeof(real_t));
    buffers[2] = malloc((local_M+2)*(local_N+2)*sizeof(real_t));

    for (int_t i = 0; i < local_M; i++)
    {
        for (int_t j = 0; j < local_N; j++)
        {
            int_t global_i = cart_coords[0] * local_M + i;
            int_t global_j = cart_coords[1] * local_N + j;
            // Calculate delta (radial distance) adjusted for M x N grid
            real_t delta = sqrt(((global_i - M/2.0) * (global_i - M/2.0)) / (real_t)M +
                                ((global_j - N/2.0) * (global_j - N/2.0)) / (real_t)N);
            U_prv(i,j) = U(i,j) = exp(-4.0*delta*delta);
        }
    }

    // Set the time step for 2D case
    dt = dx*dy / (c * sqrt(dx*dx+dy*dy));
}

// Get rid of all the memory allocations
void domain_finalize(void)
{
    free(buffers[0]);
    free(buffers[1]);
    free(buffers[2]);
}

// TASK: T5
// Integration formula
void time_step(void)
{
    for (int_t i = 0; i < local_M; i++)
    {
        for (int_t j = 0; j < local_N; j++)
        {
            U_nxt(i,j) = -U_prv(i,j) + 2.0*U(i,j)
                     + (dt*dt*c*c)/(dx*dy) * (
                        U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)-4.0*U(i,j)
                    );
        }
    }
}

// TASK: T6
// Communicate the border between processes.
void border_exchange(void)
{
    MPI_Status status;

    //MPI_Cart_shift(cart_comm, 0, 1, &north, &south);
    //MPI_Cart_shift(cart_comm, 1, 1, &west, &east);

    // Send to north, receive from south
    MPI_Sendrecv(&U(0,0), local_N, MPI_DOUBLE, north, 0,
                 &U(local_M,0), local_N, MPI_DOUBLE, south, 0,
                 cart_comm, &status);

    // Send to south, receive from north
    MPI_Sendrecv(&U(local_M-1,0), local_N, MPI_DOUBLE, south, 1,
                 &U(-1,0), local_N, MPI_DOUBLE, north, 1,
                 cart_comm, &status);

    // Send to east, receive from west
    MPI_Sendrecv(&U(0,local_N-1), local_M, MPI_DOUBLE, east, 2,
                 &U(0,-1), local_M, MPI_DOUBLE, west, 2,
                 cart_comm, &status);

    // Send to west, receive from east
    MPI_Sendrecv(&U(0,0), local_M, MPI_DOUBLE, west, 3,
                 &U(0,local_N), local_M, MPI_DOUBLE, east, 3,
                 cart_comm, &status);
}

// TASK: T7
// Neumann (reflective) boundary condition
void boundary_condition(void)
{
    if (south == MPI_PROC_NULL)
    {
        for (int_t j = 0; j < local_N; j++)
        {
            U(-1,j) = U(0,j);
        }
    }

    if (north == MPI_PROC_NULL)
    {
        for (int_t j = 0; j < local_N; j++)
        {
            U(local_M,j) = U(local_M-1,j);
        }
    }

    if (west == MPI_PROC_NULL)
    {
        for (int_t i = 0; i < local_M; i++)
        {
            U(i,-1) = U(i,0);
        }
    }

    if (east == MPI_PROC_NULL)
    {
        for (int_t i = 0; i < local_M; i++)
        {
            U(i,local_N) = U(i,local_N-1);
        }
    }
}

// TASK: T8
// Save the present time step in a numbered file under 'data/'
void domain_save(int_t step)
{
    char filename[256];
    sprintf(filename, "data/%.5ld.dat", step);
    MPI_File fh;
    MPI_File_open(cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_Offset offset = (cart_coords[0] * local_M * N + cart_coords[1] * local_N) * sizeof(real_t);
    printf("Rank %d: offset = %lld\n", world_rank, offset);
    MPI_File_set_view(fh, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

    for (int_t i = 0; i < local_M; i++)
    {
        printf("Rank %d: writing row %ld\n", world_rank, i*local_N);
        MPI_File_write(fh, &U(i,0), local_N, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&fh);
}

// Main time integration.
void simulate(void)
{
    // Go through each time step
    for (int_t iteration = 0; iteration <= max_iteration; iteration++)
    {
        if ((iteration % snapshot_freq) == 0)
        {
            domain_save(iteration / snapshot_freq);
        }

        // Derive step t+1 from steps t and t-1
        border_exchange();
        boundary_condition();
        time_step();

        // Rotate the time step buffers
        move_buffer_window();
    }
}

int main(int argc, char **argv)
{
    // TASK: T1c
    // Initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // TASK: T3
    // Distribute the user arguments to all the processes
    OPTIONS *options = NULL;
    if (world_rank == 0)
    {
        options = parse_args(argc, argv);
        if (!options)
        {
            fprintf(stderr, "Argument parsing failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        M = options->M;
        N = options->N;
        max_iteration = 500;
        snapshot_freq = options->snapshot_frequency;
    }

    MPI_Bcast(&M, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_iteration, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snapshot_freq, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

    // Create Cartesian communicator
    MPI_Dims_create(world_size, 2, cart_dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, cart_dims, cart_periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, world_rank, 2, cart_coords);

    // Set up the initial state of the domain
    domain_initialize();

    struct timeval t_start, t_end;

    // TASK: T2
    // Time your code
    MPI_Barrier(cart_comm);
    if (world_rank == 0)
    {
        gettimeofday(&t_start, NULL);
    }

    simulate();

    MPI_Barrier(cart_comm);
    if (world_rank == 0)
    {
        gettimeofday(&t_end, NULL);
        double elapsed = WALLTIME(t_end) - WALLTIME(t_start);
        printf("Elapsed time: %f seconds\n", elapsed);
    }

    // Clean up and shut down
    domain_finalize();

    // TASK: T1d
    // Finalise MPI
    MPI_Finalize();

    if (world_rank == 0 && options)
    {
        free(options);
    }

    return EXIT_SUCCESS;
}

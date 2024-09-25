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
// BEGIN: T1a
#include <mpi.h>
// END: T1a


// Convert 'struct timeval' into seconds in double prec. floating point
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

// Option to change numerical precision
typedef int64_t int_t;
typedef double real_t;


// Buffers for three time steps, indexed with 2 ghost points for the boundary
real_t
*buffers[3] = { NULL, NULL, NULL };

// TASK: T1b
// Declare variables each MPI process will need
// BEGIN: T1b
#define U_prv(i,j) buffers[0][((i)+1)*(N+2)+(j)+1]
#define U(i,j)     buffers[1][((i)+1)*(N+2)+(j)+1]
#define U_nxt(i,j) buffers[2][((i)+1)*(N+2)+(j)+1]

int world_size, world_rank;
MPI_Comm cart_comm;
int cart_dims[2], cart_period[2], cart_coords[2];

// END: T1b

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
real_t
dt;

// Add these global variables at the top of the file
int_t local_M, local_N;
int_t start_i, start_j;


// Rotate the time step buffers.
void move_buffer_window ( void )
{
    real_t *temp = buffers[0];
    buffers[0] = buffers[1];
    buffers[1] = buffers[2];
    buffers[2] = temp;
}


// TASK: T4
// Set up our three buffers, and fill two with an initial perturbation
// and set the time step.
void domain_initialize ( void )
{
    // BEGIN: T4
    int dims[2];
    MPI_Cart_get(cart_comm, 2, dims, cart_period, cart_coords);
    local_M = M / dims[0];
    local_N = N / dims[1];

    start_i = local_M * cart_coords[0];
    start_j = local_N * cart_coords[1];

    buffers[0] = malloc ( (local_M+2)*(local_N+2) * sizeof(real_t) );
    buffers[1] = malloc ( (local_M+2)*(local_N+2) * sizeof(real_t) );
    buffers[2] = malloc ( (local_M+2)*(local_N+2) * sizeof(real_t) );

    for ( int_t i=0; i<local_M; i++ )
    {
        for ( int_t j=0; j<local_N; j++ )
        {
            real_t global_i = start_i + i;
            real_t global_j = start_j + j;
            // Calculate delta (radial distance) adjusted for M x N grid
            real_t delta = sqrt(((global_i - M/2.0) * (global_i - M/2.0)) / (real_t)M +
                                ((global_j - N/2.0) * (global_j - N/2.0)) / (real_t)N);
            U_prv(i,j) = U(i,j) = exp ( -4.0*delta*delta );
        }
    }

    // Set the time step for 2D case
    dt = dx*dy / (c * sqrt (dx*dx+dy*dy));
    // END: T4
}


// Get rid of all the memory allocations
void domain_finalize ( void )
{
    free ( buffers[0] );
    free ( buffers[1] );
    free ( buffers[2] );
}


// TASK: T5
// Integration formula
void time_step ( void )
{
    // BEGIN: T5
    for ( int_t i=0; i<local_M; i++ )
    {
        for ( int_t j=0; j<local_N; j++ )
        {
            U_nxt(i,j) = -U_prv(i,j) + 2.0*U(i,j)
                + (dt*dt*c*c)/(dx*dy) * (
                U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1)-4.0*U(i,j)
            );
        }
    }
    // END: T5
}

// TASK: T6
// Communicate the border between processes.
void border_exchange ( void )
{
    // BEGIN: T6
    ;
    // END: T6
}


// TASK: T7
// Neumann (reflective) boundary condition
void boundary_condition ( void )
{
    // BEGIN: T7
    for ( int_t i=0; i<M; i++ )
    {
        U(i,-1) = U(i,1);
        U(i,N)  = U(i,N-2);
    }
    for ( int_t j=0; j<N; j++ )
    {
        U(-1,j) = U(1,j);
        U(M,j)  = U(M-2,j);
    }
    // END: T7
}


// TASK: T8
// Save the present time step in a numbered file under 'data/'
void domain_save ( int_t step )
{
    // BEGIN: T8
    char filename[256];
    sprintf ( filename, "data/%.5ld.dat", step );
    FILE *out = fopen ( filename, "wb" );
    for ( int_t i=0; i<M; i++ )
    {
        fwrite ( &U(i,0), sizeof(real_t), N, out );
    }
    fclose ( out );
    // END: T8
}


// Main time integration.
void simulate( void )
{
    // Go through each time step
    for ( int_t iteration=0; iteration<=max_iteration; iteration++ )
    {
        if ( (iteration % snapshot_freq)==0 )
        {
            domain_save ( iteration / snapshot_freq );
        }

        // Derive step t+1 from steps t and t-1
        border_exchange();
        boundary_condition();
        time_step();

        // Rotate the time step buffers
        move_buffer_window();
    }
}


int main ( int argc, char **argv )
{
    // TASK: T1c
    // Initialise MPI
    // BEGIN: T1c
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
    // END: T1c


    // TASK: T3
    // Distribute the user arguments to all the processes
    // BEGIN: T3
    OPTIONS *options = NULL;
    if (world_rank == 0) {
        options = parse_args(argc, argv);
        if (!options) {
            fprintf(stderr, "Argument parsing failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        M = options->M;
        N = options->N;
        max_iteration = options->max_iteration;
        snapshot_freq = options->snapshot_frequency;
    }

    // Broadcast the parsed options to all processes
    MPI_Bcast(&M, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_iteration, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snapshot_freq, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

    // Set up Cartesian communicator
    int dims[2] = {0, 0};  // Let MPI decide the dimensions
    int periods[2] = {0, 0};  // Non-periodic
    MPI_Dims_create(world_size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, world_rank, 2, cart_coords);
    // END: T3

    // Set up the initial state of the domain
    domain_initialize();


    struct timeval t_start, t_end;

    // TASK: T2
    // Time your code
    // BEGIN: T2
    if ( world_rank == 0 )
    {
        gettimeofday ( &t_start, NULL );
    }
    simulate();
    if ( world_rank == 0 )
    {
        gettimeofday ( &t_end, NULL );
        printf ( "Elapsed time: %.6f seconds\n", WALLTIME(t_end) - WALLTIME(t_start) );
    }
    // END: T2

    // Clean up and shut down
    domain_finalize();

    // TASK: T1d
    // Finalise MPI
    // BEGIN: T1d
    MPI_Finalize();
    // END: T1d

    exit ( EXIT_SUCCESS );
}

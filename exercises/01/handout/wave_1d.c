#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>


// Option to change numerical precision.
typedef int64_t int_t;
typedef double real_t;

// Simulation parameters: size, step count, and how often to save the state.
const int_t
    N = 1024,
    max_iteration = 4000,
    snapshot_freq = 10;

// Wave equation parameters, time step is derived from the space step.
const real_t
    c  = 1.0,
    dx  = 1.0;
real_t
    dt;

// Buffers for three time steps, indexed with 2 ghost points for the boundary.
real_t
    *buffers[3] = { NULL, NULL, NULL };


#define U_prv(i) buffers[0][(i)+1]
#define U(i)     buffers[1][(i)+1]
#define U_nxt(i) buffers[2][(i)+1]


// Save the present time step in a numbered file under 'data/'.
void domain_save ( int_t step )
{
    char filename[256];
    sprintf ( filename, "data/%.5ld.dat", step );
    FILE *out = fopen ( filename, "wb" );
    fwrite ( &U(0), sizeof(real_t), N, out );
    fclose ( out );
}


// TASK: T1
// Set up our three buffers, fill two with an initial cosine wave,
// and set the time step.
void domain_initialize ( void )
{
    // BEGIN: T1
    // Allocate memory for the buffers (N + 2 for the ghost points)
    for (int i = 0; i < 3; ++i) {
        buffers[i] = (real_t *)malloc((N + 2) * sizeof(real_t));
    }

    // Initialize the U_prv and U buffers with the initial condition
    for (int_t i = 0; i < N; ++i) {
        real_t x = (real_t)i / (real_t)(N - 1);
        U(i) = U_prv(i) = cos(M_PI * x);
    }

    // Set the time step based on the CFL condition
    dt = dx / c;
    // END: T1
}


// TASK T2:
// Return the memory to the OS.
// BEGIN: T2
void domain_finalize ( void )
{
    for (int i = 0; i < 3; ++i) {
        free(buffers[i]);
    }
}
// END: T2


// TASK: T3
// Rotate the time step buffers.
// BEGIN: T3
// TASK T3: Rotate the time step buffers.
void domain_rotate_buffers(void)
{
    real_t *temp = buffers[0];
    buffers[0] = buffers[1];
    buffers[1] = buffers[2];
    buffers[2] = temp;
}

// END: T3


// TASK: T4
// Derive step t+1 from steps t and t-1.
// BEGIN: T4
// TASK T4: Derive step t+1 from steps t and t-1.
void domain_update(void)
{
    for (int_t i = 1; i < N; i++)
    {
        U_nxt(i) = -U_prv(i) + 2 * U(i) + (dt * dt * c * c / (dx * dx)) * (U(i - 1) + U(i + 1) - 2 * U(i));
    }
}

// END: T4


// TASK: T5
// Neumann (reflective) boundary condition.
// BEGIN: T5
// TASK T5: Neumann (reflective) boundary condition.
void domain_apply_boundary(void)
{
    // Reflect the value at the boundary by mirroring the interior value.
    U(0) = U(1);       // Left boundary
    U(N) = U(N-1);   // Right boundary
}

// END: T5


// TASK: T6
// Main time integration.
void simulate( void )
{
// BEGIN: T6
    for (int_t iteration = 0; iteration < max_iteration; iteration++)
    {
        // Apply Neumann boundary conditions
        domain_apply_boundary();

        // Update the wave equation for the next time step
        domain_update();

        // Rotate the buffers to advance in time
        domain_rotate_buffers();

        // Save the state at specified intervals
        if (iteration % snapshot_freq == 0)
        {
            domain_save(iteration/snapshot_freq);
        }
    }
// END: T6
}


int main ( void )
{
    domain_initialize();

    simulate();

    domain_finalize();
    exit ( EXIT_SUCCESS );
}

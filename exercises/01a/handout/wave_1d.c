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
    // Allocate memory for three sets of spatial domains: the previous, the current, and the next domain.
    // The domain includes two ghost points at each end.
    buffers[0] = (real_t *)malloc((N + 2) * sizeof(real_t));
    buffers[1] = (real_t *)malloc((N + 2) * sizeof(real_t));
    buffers[2] = (real_t *)malloc((N + 2) * sizeof(real_t));

    // Apply the initial condition from Equation 4 to the current and previous time steps.
    for (int_t i = 0; i < N; ++i) {
        real_t x = (real_t)i / (real_t)(N - 1);
        U(i) = U_prv(i) = cos(M_PI * x);
    }

    // Set the time step, dt, based on Equation 5.
    dt = dx / c;
    // END: T1
}


// TASK T2:
// Return the memory to the OS.
// BEGIN: T2
void domain_finalize ( void )
{
    // Implement the function domain finalize. The function should free the memory allocated by domain initialize.
    free(buffers[0]);
    free(buffers[1]);
    free(buffers[2]);
}
// END: T2


// TASK: T3
// Rotate the time step buffers.
// BEGIN: T3
void domain_rotate_buffers(void)
{
    // Define and implement a function that shifts the buffer domain window one step forward.
    real_t *temp = buffers[0]; // Save the previous buffer
    buffers[0] = buffers[1]; // Move the current buffer to the previous buffer
    buffers[1] = buffers[2]; // Move the next buffer to the current buffer
    buffers[2] = temp; // Move the previous buffer to the next buffer so that it can be reused
}
// END: T3


// TASK: T4
// Derive step t+1 from steps t and t-1.
// BEGIN: T4
void domain_update(void)
{
    // Define and implement a function that does one iteration of integration over the time domain (does one step forward in time).
    // The function should calculate the next time iteration based on Equation 6 for every point in space.
    // Used the equation shown in the slides of recitation 2 to update the wave equation.
    for (int_t i = 1; i < N; i++)
    {
        U_nxt(i) = -U_prv(i) + 2 * U(i) + (dt * dt * c * c / (dx * dx)) * (U(i - 1) + U(i + 1) - 2 * U(i));
    }
}
// END: T4


// TASK: T5
// Neumann (reflective) boundary condition.
// BEGIN: T5
void domain_apply_boundary(void)
{
    // Define and implement a function that sets the value of the ghost cell to mirror that of the point one further into the domain.
    // Reflect the value at the boundary by mirroring the interior value.
    // From my understanding, to implement the Neumann boundary condition, we need to have the equation U(0) = U(1) and U(N) = U(N-1), where:
    // U(0) is the value at the left boundary (ghost point)
    // U(1) is the value at the first point in the domain
    // U(N) is the value at the right boundary (ghost point)
    // U(N-1) is the value at the last point in the domain
    U(0) = U(1);
    U(N) = U(N-1);
}
// END: T5


// TASK: T6
// Main time integration.
void simulate( void )
{
// BEGIN: T6
    for (int_t iteration = 0; iteration < max_iteration; iteration++)
    {
        // Apply the boundary condition to the ghost cells (T5).
        domain_apply_boundary();

        //  Do one time step forward (T4).
        domain_update();
        domain_rotate_buffers();

        // Store the domain in a .dat file if the time iteration matches that of the snapshot frequency, snapshot freq, by using the domain save function.
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

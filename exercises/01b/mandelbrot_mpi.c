#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define XSIZE 2560
#define YSIZE 2048

#define MAXITER 512 

double xleft = -2.01;
double xright = 1;
double yupper, ylower;
double ycenter = 1e-6;
double step;

int pixel[XSIZE * YSIZE];

#define PIXEL(i, j) ((i) + (j) * XSIZE)

typedef struct {
    double real, imag;
} complex_t;

void calculate(int start_row, int end_row, int rank) {
    printf("Process %d: Starting calculation for rows %d to %d\n", rank, start_row, end_row - 1);
    
    for (int i = 0; i < XSIZE; i++) {
        for (int j = start_row; j < end_row; j++) {
            /* Calculate the number of iterations until divergence for each pixel.
               If divergence never happens, return MAXITER */
            complex_t c, z, temp;
            int iter = 0;
            c.real = (xleft + step * i);
            c.imag = (ylower + step * j);
            z = c;

            while (z.real * z.real + z.imag * z.imag < 4) {
                temp.real = z.real * z.real - z.imag * z.imag + c.real;
                temp.imag = 2 * z.real * z.imag + c.imag;
                z = temp;
                if (++iter == MAXITER) break;
            }

            pixel[PIXEL(i, j)] = iter;
        }
        
    }
    
    printf("Process %d: Finished calculation for rows %d to %d\n", rank, start_row, end_row - 1);
}

typedef unsigned char uchar;

// Function to save the image as a BMP file
void savebmp(char *name, uchar *buffer, int x, int y) {
    FILE *f = fopen(name, "wb");
    if (!f) {
        printf("Error writing image to disk.\n");
        return;
    }
    unsigned int size = x * y * 3 + 54;
    uchar header[54] = {'B', 'M', size & 255, (size >> 8) & 255, (size >> 16) & 255, size >> 24, 0,
                                0, 0, 0, 54, 0, 0, 0, 40, 0, 0, 0, x & 255, x >> 8, 0, 0, y & 255, y >> 8, 0, 0, 1, 0, 24, 0,
                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    fwrite(header, 1, 54, f);
    fwrite(buffer, 1, XSIZE * YSIZE * 3, f);
    fclose(f);
}

// Function to set color based on iteration count
void fancycolour(unsigned char *p, int iter) {
    if (iter == MAXITER);
    else if (iter < 8) { p[0] = 128 + iter * 16; p[1] = p[2] = 0; }
    else if (iter < 24) { p[0] = 255; p[1] = p[2] = (iter - 8) * 16; }
    else if (iter < 160) { p[0] = p[1] = 255 - (iter - 24) * 2; p[2] = 255; }
    else { p[0] = p[1] = (iter - 160) * 2; p[2] = 255 - (iter - 160) * 2; }
}


int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 3) {
        if (rank == 0) {
            printf("This program should be run with 3 processes.\n");
        }
        MPI_Finalize();
        return 1;
    }

    /* Calculate the range in the y-axis such that we preserve the
       aspect ratio */
    step = (xright - xleft) / XSIZE;
    yupper = ycenter + (step * YSIZE) / 2;
    ylower = ycenter - (step * YSIZE) / 2;

    // Determine the rows each process will handle
    int rows_per_process = YSIZE / size;
    int start_row = rank * rows_per_process;
    int end_row = (rank == size - 1) ? YSIZE : start_row + rows_per_process;

    printf("Process %d: Assigned rows %d to %d\n", rank, start_row, end_row - 1);

    double start_time = MPI_Wtime();

    // Perform the calculation for the assigned rows
    calculate(start_row, end_row, rank);

    if (rank != 0) {
        printf("Process %d: Sending data to Process 0\n", rank);
        MPI_Send(&pixel[PIXEL(0, start_row)], (end_row - start_row) * XSIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
        printf("Process %d: Data sent to Process 0\n", rank);
    } else {
        // Process 0 receives data from other processes
        for (int i = 1; i < size; i++) {
            int recv_start_row = i * rows_per_process;
            int recv_end_row = (i == size - 1) ? YSIZE : recv_start_row + rows_per_process;
            printf("Process 0: Receiving data from Process %d\n", i);
            MPI_Recv(&pixel[PIXEL(0, recv_start_row)], (recv_end_row - recv_start_row) * XSIZE, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process 0: Received data from Process %d\n", i);
        }

        double end_time = MPI_Wtime();
        printf("Total calculation time: %f seconds\n", end_time - start_time);

        printf("Process 0: Creating image buffer\n");
        unsigned char *buffer = calloc(XSIZE * YSIZE * 3, 1);
        for (int i = 0; i < XSIZE; i++) {
            for (int j = 0; j < YSIZE; j++) {
                int p = ((YSIZE - j - 1) * XSIZE + i) * 3;
                fancycolour(buffer + p, pixel[PIXEL(i, j)]);
            }
        }

        printf("Process 0: Saving image to disk\n");
        savebmp("mandel_mpi.bmp", buffer, XSIZE, YSIZE);
        free(buffer);
        printf("Process 0: Image saved as 'mandel_mpi.bmp'\n");
    }

    MPI_Finalize();
    return 0;
}

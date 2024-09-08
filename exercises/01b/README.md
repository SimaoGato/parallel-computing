# Exercise 1b - Mandelbrot
This is a quick README for the Mandelbrot exercise.

## Structure
This repository contains the following:
- 'mandelbrot_mpi.c': The C code for the Mandelbrot exercise, parallelized using MPI.
- 'Makefile': The Makefile for the mandelbrot_mpi.c code.
- '01b-answers.pdf': The answers to the theory questions in the exercise.
- 'mandel-serial/': Directory containing the given serial version of the Mandelbrot code and the Makefile for it.

## Compile and Run the MPI Code
To compile the code, run the following command:
```bash
make
```

To run the code, use the following command:
```bash
make run
```
This command will already run the code with 3 processes, as requested in the exercise. It produces a image file called 'mandel_mpi.bmp'.

To clean the directory, run the following command:
```bash
make clean
```

# Compiler and flags
MPICC = mpicc
CFLAGS = -Wall -O3
LDFLAGS = -lm

# Executable name
TARGET = mandelbrot_mpi

# Source file
SRC = mandelbrot_mpi.c

# Number of processes to run
NP = 3

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SRC)
	$(MPICC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Run rule
run: $(TARGET)
	mpirun -np $(NP) ./$(TARGET)

# Clean up
clean:
	rm -f $(TARGET) mandel_mpi.bmp

# Phony targets
.PHONY: all run clean

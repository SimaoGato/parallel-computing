# Compiler and flags
CC = gcc
CFLAGS = -O2 -Wall -std=c99 -lm

# Executable name
TARGET = mandelbrot

# Source file
SRC = mandelbrot.c

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

# Run rule
run: $(TARGET)
	./$(TARGET) 1

# Clean up
clean:
	rm -f $(TARGET) mandel.bmp

# Phony targets
.PHONY: all run clean

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -O3

# Target executable
TARGET = myprogram

# Source files
SRC = main.c bitmap.c

# Output files to be cleaned
CLEAN_FILES = $(TARGET) after.bmp

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

# Clean up generated files
clean:
	rm -f $(CLEAN_FILES)

# Phony targets (they aren't files)
.PHONY: all clean

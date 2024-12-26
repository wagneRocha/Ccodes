# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -std=c11 -O3 -march=native -flto=auto -g
# CFLAGS = -Wall -Wextra -std=c11 -O0 -flto -g

# Linker flags
LDFLAGS = -llapacke -llapack -lblas -lm -flto

# Target executable
TARGET = main

# Source files
SRCS = ADT_matrices.c ADT_vectors.c ADT_stack.c ADT_strings.c ADT_files.c ADT_DDGP.c main.c

# Object files
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJS)

# Phony targets
.PHONY: all clean

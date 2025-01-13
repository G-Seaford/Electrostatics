# Makefile for compiling netCDF_Interface module and test program

# Compiler
FC = gfortran

# List of Fortran source files
PROGRAM_FILES = common_data_structure.f90 command_line.f90 error_handler.f90 velocity_verlet.f90 netcdffix_array.f90 gauss_seidel.f90 Main.f90

# Name of the compiled executable
OUTFILE = test

# Compiler flags: Use nf-config to grab the compile flags
FFLAGS = -g $(shell nf-config --fflags) -std=f2018 -Wall -fcheck=all

# Linker flags: Use nf-config to grab the link flags
FLIBS = $(shell nf-config --flibs)

# Detect Operating System
UNAME_S := $(shell uname -s)

# Remove duplicate '-lnetcdf' on macOS if nf-config provides it twice
ifeq ($(UNAME_S),Darwin)
    # Remove the second occurrence of '-lnetcdf'
    FLIBS := $(firstword $(FLIBS)) $(filter-out -lnetcdf,$(FLIBS))
endif

# Derived object files
OBJECTS = $(PROGRAM_FILES:.f90=.o)

# Default target
all: $(OUTFILE)

# Rule to build the executable
$(OUTFILE): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) $(FLIBS) -o $(OUTFILE)

# Generic rule for compiling .f90 files to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up compiled objects and executable
clean:
	rm -f $(OBJECTS) $(OUTFILE)

# Phony targets
.PHONY: all clean


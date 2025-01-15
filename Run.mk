# Makefile for compiling the Gauss-Seidel data with command line arguments,
# running the visualisation code before cleaning up executables.
# @author Gianluca Seaford

#####################################################
# 				Simulation Parameters				#
#####################################################

command_line_arguments := nx=100 ny=100 problem=double
#####################################################

# Compiler
FC = gfortran

# List of Fortran source files
PROGRAM_FILES = common_data_structure.f90 command_line.f90 error_handler.f90 velocity_verlet.f90 netcdffix_array.f90 gauss_seidel.f90 Main.f90

# Name of the compiled executable
OUTFILE = Electrostatics

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

.PHONY: default
default: RunAllFiles

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

##################################################################
## 				Script to run all files in order 				##
##################################################################

.PHONY: RunAllFiles
RunAllFiles: all
	@echo ">>> 	Running the executable with command-line arguments... 	<<<"
	./$(OUTFILE) $(command_line_arguments)
	@echo ">>>															<<<"

	@echo ">>> 				Calling the Python script...	 			<<<"
	python3 vis.py
	@echo ">>>															<<<"

	@echo ">>> 					Cleaning up artifacts 					<<<"
	rm -f $(OBJECTS) $(OUTFILE)
	@echo ">>>															<<<"
	@echo ">>> 				Process terminated successfully 			<<<"
	@echo ">>>															<<<"

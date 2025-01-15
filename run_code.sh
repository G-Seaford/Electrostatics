#!/bin/bash

# Bash script to implement the workflow for the Gauss-Seidel based electrostatics problem.
# @author Gianluca Seaford


# Exit immediately on any error
set -e

#########################
#   Default Arguments   #
#########################
nx=100
ny=100
problem=single

extra_args=()
help_requested=0

# Allow the passing of command line arguments with the bash script.
for arg in "$@"; do
  case $arg in
    nx=*)
      nx="${arg#*=}"
      ;;
    ny=*)
      ny="${arg#*=}"
      ;;
    problem=*)
      problem="${arg#*=}"
      ;;
    --help)
      help_requested=1
      extra_args+=("$arg")  # Allows the identification of the '--help' command line flag.
      ;;
    *)
      extra_args+=("$arg") # Pass extra arguments to Fortran for handling.
      ;;
  esac
done


echo ">>> 					Compiling Fortran Code					<<<"
make  -f compile_fortran.mk 

echo ">>> 	Running the executable with command-line arguments... 	<<<"

if [ "$help_requested" = "1" ]; then
    echo ">>>     --help requested. Displaying help then exiting.       <<<"
    ./Electrostatics --help
    echo ">>>															<<<"
    echo ">>> 					Cleaning up artifacts 					<<<"
    make clean -f compile_fortran.mk
    exit 0
fi


./Electrostatics nx="${nx}" ny="${ny}" problem="${problem}" "${extra_args[@]}"
echo ">>>															<<<"

echo ">>> 				Calling the Python script        			<<<"
python3 vis.py
echo ">>>															<<<"

echo ">>> 					Cleaning up artifacts 					<<<"
make clean -f compile_fortran.mk
echo ">>>															<<<"

echo ">>> 				Process terminated successfully 			<<<"
echo ">>>															<<<"

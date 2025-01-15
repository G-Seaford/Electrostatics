# Electrostatics

**Electrostatics Repository for PX913**  
**Authors**: Gianluca Seaford & Facundo Costa

This project implements a Gauss–Seidel solver for electrostatic problems in Fortran, along with a Python script for visualising the results. A makefile is used to compile the Fortran code, and a bash script (`run_code.sh`) is used to implement the workflow.

---

## Requirements

- **Fortran compiler** (e.g., `gfortran`) **with netCDF support**  
- **Python 3.x** with required libraries (e.g., `matplotlib`, `numpy`, `netCDF4`)  
- **Make** utility (to build the Fortran code)

---

## Files

- **`compile_fortran.mk`**  
  A Makefile that defines how to compile and link the Fortran source files into an executable called `Electrostatics`. It also provides a `clean` target to remove build artefacts.

- **`run_code.sh`**  
  A bash script that:
  1. Compiles the Fortran code by calling:
     ```bash
     make -f compile_fortran.mk
     ```
  2. Runs the compiled executable with various command-line arguments (including defaults for `nx=100`, `ny=100`, `problem=single`).
  3. Invokes:
     ```bash
     python3 vis.py
     ```
     for visualisation.
  4. Cleans up build artefacts by calling:
     ```bash
     make clean -f compile_fortran.mk
     ```

- **`vis.py`**  
  A Python script for visualising the output (e.g., reading `.nc` files and producing plots).

---

## Usage
** Compile & Run with Defaults: **
        ./run_code.sh

    This will compile the code and run './Electrostatics nx=100 ny=100 problem=single' before calling the python3 'vis.py' visualisation script and cleaning up any artifacts from the executable construction.

    * Request Help:
        ./run_code.sh --help

    For help with command line arguments run the above line. The script will compile the code, run './Electrostatics --help' before displaying the help window and cleaning up any artifacts from the executable construction.

    * Custom Parameters:

    You can override any or all of the default parameters, and add additional parameters. The parameters available are given by the '--help' menu.

    * Cleaning Up:

    In normal operation, the script automatically cleans when it completes. If you simply wish to remove compiled files manually, you can run 'make -f compile_fortran.mk clean' to remove the Electrostatics executable and any artifacts from the executable construction.
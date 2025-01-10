MODULE run_data_module
  ! This module contains the run data that I will then put on the file
  ! which store all the important details of the simulation that was run

  use ISO_FORTRAN_ENV

  IMPLICIT NONE

  TYPE, PUBLIC :: run_data
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: axis_x
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: axis_y
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: E_x
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: E_y
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: position
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: velocity
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: acceleration
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: potential
    INTEGER :: num_steps
    INTEGER :: nx
    INTEGER :: ny
  END TYPE run_data

END MODULE run_data_module

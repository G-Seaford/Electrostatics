MODULE CommonDataStructure
    !> Module for defining the common data structure.
    !! @author Gianluca Seaford
    !! @date 2021-06-01
    !! @version 1.0

    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    TYPE :: CommonData
        !> Derived type for data storage across modules

        !! Simulation Parameters
        INTEGER         :: nx, ny, n_iter
        REAL(REAL64)    :: dx, dy, dt

        !! Poisson equation data
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE   :: phi
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE   :: rho

        !! Electric field data
        REAL(REAL64), DIMENSION(:, :), ALLOCATABLE  :: Ex
        REAL(REAL64), DIMENSION(:, :), ALLOCATABLE  :: Ey

        !! Electron data
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE     :: positions
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE     :: velocities
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE     :: accelerations

    END TYPE CommonData

END MODULE CommonDataStructure
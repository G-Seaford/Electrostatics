MODULE solve_gauss_seidel
    !> Solves the Poisson equation iteratively via the Gauss-Seidel method.
    !! author: Gianluca Seaford
    !! version: 1.0
    !! date: 2024-11-29
    
    IMPLICIT NONE

    !> Define derived type for data storage
    TYPE :: GaussSeidelData
        INTEGER :: nx, ny
        REAL :: dx, dy

        !! Poisson equation data
        REAL, DIMENSION(:,:), ALLOCATABLE :: phi
        REAL, DIMENSION(:,:), ALLOCATABLE :: rho

        !! Electric field data
        REAL, DIMENSION(:), ALLOCATABLE :: Ex
        REAL, DIMENSION(:), ALLOCATABLE :: Ey

        !! Electron data
        REAL, DIMENSION(:), ALLOCATABLE :: position_x, position_y
        REAL, DIMENSION(:), ALLOCATABLE :: velocity_x, velocity_y
        REAL, DIMENSION(:), ALLOCATABLE :: acceleration_x, acceleration_y

    END TYPE GaussSeidelData

    CONTAINS


END MODULE solve_gauss_seidel
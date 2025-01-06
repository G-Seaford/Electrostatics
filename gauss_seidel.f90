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

    SUBROUTINE InitialiseGaussSeidelData(data, init_type, nx, ny, dx, dy, status)
        !> Initialises the data structure for the Gauss-Seidel method.
        !! @param data: The GaussSeidelData structure to initialise.
        !! @param nx: The number of grid points in the x-direction.
        !! @param ny: The number of grid points in the y-direction.
        !! @param dx: The grid spacing in the x-direction.
        !! @param dy: The grid spacing in the y-direction.

        TYPE(GaussSeidelData), INTENT(INOUT) :: data
        CHARACTER(LEN=*), INTENT(IN) :: init_type
        INTEGER, INTENT(IN) :: nx, ny
        REAL, INTENT(IN) :: dx, dy
        INTEGER, INTENT(INOUT) :: status

        data%nx = nx
        data%ny = ny

        data%dx = dx
        data%dy = dy

        ALLOCATE(data%phi(0:nx+1, 0:ny+1), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for phi"
            RETURN
        END IF

        ALLOCATE(data%rho(0:nx+1, 0:ny+1), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for rho"
            RETURN
        END IF

        SELECT CASE(init_type)

            CASE("null")
                !! Initialise phi and rho to zero
                data%phi = 0.0
                data%rho = 0.0

            CASE("single")

            CASE("double")
        
        END SELECT

        

    END SUBROUTINE InitialiseGaussSeidelData

END MODULE solve_gauss_seidel
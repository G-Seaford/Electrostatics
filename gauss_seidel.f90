MODULE solve_gauss_seidel
    !> Solves the Poisson equation iteratively via the Gauss-Seidel method.
    !! @author Gianluca Seaford
    !! @version 1.0
    !! @date 2024-11-29
    
    IMPLICIT NONE

    TYPE :: GaussSeidelData
        !> Derived type for data storage across modules

        !! Simulation Parameters
        INTEGER         :: nx, ny, n_iter
        REAL            :: dx, dy, dt

        !! Poisson equation data
        REAL, DIMENSION(:,:), ALLOCATABLE   :: phi
        REAL, DIMENSION(:,:), ALLOCATABLE   :: rho

        !! Electric field data
        REAL, DIMENSION(:, :), ALLOCATABLE  :: Ex
        REAL, DIMENSION(:, :), ALLOCATABLE  :: Ey

        !! Electron data
        REAL, DIMENSION(:), ALLOCATABLE     :: position_x, position_y
        REAL, DIMENSION(:), ALLOCATABLE     :: velocity_x, velocity_y
        REAL, DIMENSION(:), ALLOCATABLE     :: acceleration_x, acceleration_y

    END TYPE GaussSeidelData

    CONTAINS

    SUBROUTINE InitialiseGaussSeidelData(data, init_type, nx, ny, n_iter, dx, dy, dt, status)
        !> Initialises the data structure for the Gauss-Seidel method.
        !! @param[inout] data The GaussSeidelData structure to initialise.
        !! @param[in] init_type The initialisation type.
        !! @param[in] nx The number of grid points in the x-direction.
        !! @param[in] ny The number of grid points in the y-direction.
        !! @param[in] n_iter The number of iterations.
        !! @param[in] dx The grid spacing in the x-direction.
        !! @param[in] dy The grid spacing in the y-direction.
        !! @param[in] dt The time step.
        !! @param[inout] status Status for error handling.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data
        CHARACTER(LEN=*), INTENT(IN)            :: init_type
        INTEGER, INTENT(IN)                     :: nx, ny, n_iter
        INTEGER                                 :: i, j             !! Loop variables
        REAL, INTENT(IN)                        :: dx, dy, dt
        REAL                                    :: x, y             !! Position variables
        INTEGER, INTENT(INOUT)                  :: status

        data%nx = nx
        data%ny = ny
        data%n_iter = n_iter

        data%dx = dx
        data%dy = dy
        data%dt = dt

        !> Allocate field arrays
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

        ALLOCATE(data%Ex(0:nx+1, 0:ny+1), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for Ex"
            RETURN
        END IF

        ALLOCATE(data%Ey(0:nx+1, 0:ny+1), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for Ey"
            RETURN
        END IF

        !> Allocate particle storage arrays.
        ALLOCATE(data%position_x(0:n_iter), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for x positions"
            RETURN
        END IF

        ALLOCATE(data%position_y(0:n_iter), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for y positions"
            RETURN
        END IF

        ALLOCATE(data%velocity_x(0:n_iter), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for x velocities"
            RETURN
        END IF
        
        ALLOCATE(data%velocity_y(0:n_iter), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for y velocities"
            RETURN
        END IF

        ALLOCATE(data%acceleration_x(0:n_iter), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for x accelerations"
            RETURN
        END IF
        
        ALLOCATE(data%acceleration_y(0:n_iter), STAT=status)
        IF (status /= 0) THEN
            PRINT *, "Error allocating memory for y accelerations"
            RETURN
        END IF

        !> Initialise arrays with zero values
        data%phi = 0.0
        data%rho = 0.0

        data%Ex = 0.0
        data%Ey = 0.0

        data%position_x = 0.0
        data%position_y = 0.0
        data%velocity_x = 0.0
        data%velocity_y = 0.0
        data%acceleration_x = 0.0
        data%acceleration_y = 0.0

        SELECT CASE(init_type)

            CASE("null")
                !> Zero Gaussian peaks, with null rho

                data%velocity_x(0) = 0.1
                data%velocity_y(0) = 0.1

            CASE("single")
                !> Single Gaussian Peak

                DO j = 1, data%ny
                    y = -1.0 + (j - 1) * 2.0 / (data%ny - 1)

                    DO i = 1, data%nx
                        x = -1.0 + (i - 1) * 2.0 / (data%nx - 1)
        
                        data%rho(i, j) = EXP( - ( (x/0.1)**2 + (y/0.1)**2 ) )
                    END DO
                END DO        

                data%position_x(0) = 0.1

            CASE("double")
                !> Double Gaussian peak

                DO j = 1, data%ny
                    y = -1.0 + (j - 1) * 2.0 / (data%ny - 1)

                    DO i = 1, data%nx
                        x = -1.0 + (i - 1) * 2.0 / (data%nx - 1)
        
                        data%rho(i, j) = EXP(-(((x+0.25)/0.1)**2 + ((y+0.25)/0.1)**2))  &
                                    &  + EXP(-(((x-0.75)/0.2)**2 + ((y-0.75)/0.2)**2))
                    END DO
                END DO
        
                data%position_y = 0.5

        END SELECT

    END SUBROUTINE InitialiseGaussSeidelData

    SUBROUTINE PhiUpdate(data_in, x_idx, y_idx)
        !> Subroutine to update a given phi element
        !! @param[inout] data_in The GaussSeidelData structure holding data.
        !! @param[in] x_idx The x index.
        !! @param[in] y_idx The y index.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER, INTENT(IN)                     :: x_idx, y_idx
        REAL                                    :: dx, dy, dx2_inv, dy2_inv, rho, phi_current, phi_xdiff, phi_ydiff, denom_inv


        dx = data_in%dx
        dy = data_in%dy
        rho = data_in%rho(x_idx, y_idx)
        phi_current = data_in%phi(x_idx, y_idx)

        dx2_inv = 1/(dx * dx)
        dy2_inv = 1/(dy * dy)
        denom_inv = 1/(2 * (dx2_inv + dy2_inv))

        phi_xdiff = (data_in%phi(x_idx + 1, y_idx) + data_in%phi(x_idx - 1, y_idx)) * dx2_inv
        phi_ydiff = (data_in%phi(x_idx, y_idx + 1) + data_in%phi(x_idx, y_idx - 1)) * dy2_inv

        data_in%phi(x_idx, y_idx) = -1.0 * denom_inv * (rho - phi_xdiff - phi_ydiff)

    END SUBROUTINE PhiUpdate

    SUBROUTINE EfieldUpdate(data_in, x_idx, y_idx)
        !> Updates the electric field components at a given grid location.
        !! @param[inout] data_in The GaussSeidelData structure holding data.
        !! @param[in] x_idx The x index.
        !! @param[in] y_idx The y index.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER, INTENT(IN)                     :: x_idx, y_idx
        REAL                                    :: dx, dy, phi_xdiff, phi_ydiff

        dx = data_in%dx
        dy = data_in%dy

        phi_xdiff = (data_in%phi(x_idx + 1, y_idx) - data_in%phi(x_idx - 1, y_idx))
        phi_ydiff = (data_in%phi(x_idx, y_idx + 1) - data_in%phi(x_idx, y_idx - 1))

        data_in%Ex(x_idx, y_idx) = phi_xdiff / (2 * dx)
        data_in%Ey(x_idx, y_idx) = phi_ydiff / (2 * dy)

    END SUBROUTINE EfieldUpdate

    SUBROUTINE SweepPhi(data_in)
        !> Performs a Gauss–Seidel sweep over the entire phi array,
        !! updating phi and the electric field at each grid point.
        !! @param[inout] data_in The GaussSeidelData structure containing the grid and field arrays.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER                                 :: nx, ny, i_idx, j_idx

        DO j_idx=1, ny
            DO i_idx=1, nx

                CALL PhiUpdate(data_in, i_idx, j_idx)
                CALL EfieldUpdate(data_in, i_idx, j_idx)

            END DO

        END DO

    END SUBROUTINE SweepPhi

END MODULE solve_gauss_seidel
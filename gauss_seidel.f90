MODULE solve_gauss_seidel
    !> Solves the Poisson equation iteratively via the Gauss-Seidel method.
    !! @author Gianluca Seaford
    !! @version 1.0
    !! @date 2024-11-29

    USE command_line
    USE Error_Logging
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Insert used packages !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!

    TYPE :: GaussSeidelData
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
        REAL(REAL64), DIMENSION(:), ALLOCATABLE     :: position_x, position_y
        REAL(REAL64), DIMENSION(:), ALLOCATABLE     :: velocity_x, velocity_y
        REAL(REAL64), DIMENSION(:), ALLOCATABLE     :: acceleration_x, acceleration_y

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
        REAL(REAL64), INTENT(IN)                :: dx, dy, dt
        REAL(REAL64)                            :: x, y             !! Position variables
        INTEGER, INTENT(INOUT)                  :: status

        !> Set parameters within GaussSeidelData structure.
        data%nx = nx
        data%ny = ny
        data%n_iter = n_iter

        data%dx = dx
        data%dy = dy
        data%dt = dt

        !> Allocate field arrays.
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

        !> Initialise arrays with zero values.
        !! This is equivalent to imposing Dirichlet conditions at boundaries.
        data%phi = 0.0_REAL64
        data%rho = 0.0_REAL64

        data%Ex = 0.0_REAL64
        data%Ey = 0.0_REAL64

        data%position_x = 0.0_REAL64
        data%position_y = 0.0_REAL64
        data%velocity_x = 0.0_REAL64
        data%velocity_y = 0.0_REAL64
        data%acceleration_x = 0.0_REAL64
        data%acceleration_y = 0.0_REAL64

        SELECT CASE(init_type)

            CASE("null")
                !> Zero Gaussian peaks, with null rho.
                !! Used for a lone charge.

                data%velocity_x(0) = 0.1_REAL64
                data%velocity_y(0) = 0.1_REAL64

            CASE("single")
                !> Single Gaussian Peak.
                !! Used for a localised point charge distribution.

                DO j = 1, data%ny
                    y = -1.0 + (j - 1) * 2.0 / (data%ny - 1)

                    DO i = 1, data%nx
                        x = -1.0 + (i - 1) * 2.0 / (data%nx - 1)
        
                        data%rho(i, j) = EXP( - ( (x/0.1)**2 + (y/0.1)**2 ) )
                    END DO
                END DO        

                data%position_x(0) = 0.1_REAL64

            CASE("double")
                !> Double Gaussian peak.
                !! Used for two localised point charge distributions. 

                DO j = 1, data%ny
                    y = -1.0 + (j - 1) * 2.0 / (data%ny - 1)

                    DO i = 1, data%nx
                        x = -1.0 + (i - 1) * 2.0 / (data%nx - 1)
        
                        data%rho(i, j) = EXP(-(((x+0.25)/0.1)**2 + ((y+0.25)/0.1)**2))  &
                                    &  + EXP(-(((x-0.75)/0.2)**2 + ((y-0.75)/0.2)**2))
                    END DO
                END DO
        
                data%position_y = 0.5_REAL64

        END SELECT

    END SUBROUTINE InitialiseGaussSeidelData

    SUBROUTINE PhiUpdate(data_in, x_idx, y_idx)
        !> Subroutine to update a given phi element.
        !! Calculates the second-order partial derivatives via a central finite difference method.
        !! @param[inout] data_in The GaussSeidelData structure holding data.
        !! @param[in] x_idx The x index.
        !! @param[in] y_idx The y index.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER, INTENT(IN)                     :: x_idx, y_idx
        REAL(REAL64)                            :: dx, dy, dx2_inv, dy2_inv, rho, phi_current, phi_xdiff, phi_ydiff, denom_inv


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

    SUBROUTINE EfieldUpdate(data_in)
        !> Subroutine to update the electric field components in Ex and Ey.
        !! Uses a modified first-order 
        !! @param[inout] data_in The GaussSeidelData structure holding data.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER                                 :: x_idx, y_idx
        REAL(REAL64)                            :: dx, dy, phi_xdiff, phi_ydiff

        dx = data_in%dx
        dy = data_in%dy

        DO y_idx = 1, data_in%ny
            DO x_idx = 1, data_in%nx

                phi_xdiff = (data_in%phi(x_idx + 1, y_idx) - data_in%phi(x_idx - 1, y_idx))
                phi_ydiff = (data_in%phi(x_idx, y_idx + 1) - data_in%phi(x_idx, y_idx - 1))

                data_in%Ex(x_idx, y_idx) = phi_xdiff / (2 * dx)
                data_in%Ey(x_idx, y_idx) = phi_ydiff / (2 * dy)

            END DO
        END DO

    END SUBROUTINE EfieldUpdate

    LOGICAL FUNCTION ConvergenceTest(data_in) RESULT(is_converged)
        !> Function to test whether the convergence criterion is achieved.
        !! Uses the ratio of the total error to rms distance to evaluate convergence.
        !! @param[in] data_in The GaussSeidelData structure holding data.
        !! @param[out] isconverged Logical determining whether convergence has been attained.

        TYPE(GaussSeidelData), INTENT(IN)   :: data_in
        INTEGER                             :: x_idx, y_idx, n_sites
        REAL(REAL64)                        :: rho, dx, dy, dx2_inv, dy2_inv, phi_xdiff, phi_ydiff, tot_diff, e_tot, d_rms, ratio, epsilon = 1.0e-12_REAL64

        is_converged = .FALSE.

        dx = data_in%dx
        dy = data_in%dy
        n_sites = data_in%nx * data_in%ny

        e_tot = 0.0_REAL64
        d_rms = 0.0_REAL64

        dx2_inv = 1/(dx * dx)
        dy2_inv = 1/(dy * dy)

        DO y_idx = 1, data_in%ny
            DO x_idx = 1, data_in%nx
                !> Loop over (x,y) pairs and calculate the rms square distance and total error.

                rho = data_in%rho(x_idx, y_idx)

                phi_xdiff = (data_in%phi(x_idx + 1, y_idx) - 2.0 * data_in%phi(x_idx, y_idx) + data_in%phi(x_idx - 1, y_idx))
                phi_ydiff = (data_in%phi(x_idx, y_idx + 1) -  2.0 * data_in%phi(x_idx, y_idx) + data_in%phi(x_idx, y_idx - 1))

                tot_diff = phi_xdiff * dx2_inv + phi_ydiff * dy2_inv

                d_rms = d_rms + tot_diff
                e_tot = e_tot + ABS(tot_diff - rho)

            END DO
        END DO

        !> Normalise the rms squared distance and sqrt.
        d_rms = SQRT(d_rms / n_sites)

        IF (d_rms < epsilon) THEN
            !> Prevent division by zero.
            IF (e_tot < epsilon) THEN
                !> True convergence achieved.
                ratio = 0.0_REAL64
            
            ELSE
                !> Indicative of an ill-conditioned Laplacian matrix
                ratio = HUGE(1.0_REAL64)
            END IF 

        ELSE 
            ratio = e_tot / d_rms

        END IF

        IF (ratio < 0.00001_REAL64) THEN
            is_converged = .TRUE.
        END IF

    END FUNCTION ConvergenceTest

    SUBROUTINE SweepPhi(data_in)
        !> Performs a Gauss–Seidel sweep over the entire phi array,
        !! updating phi at each grid point.
        !! @param[inout] data_in The GaussSeidelData structure containing the grid and field arrays.

        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER                                 :: nx, ny, x_idx, y_idx

        DO y_idx=1, ny
            DO x_idx=1, nx

                CALL PhiUpdate(data_in, x_idx, y_idx)

            END DO

        END DO

    END SUBROUTINE SweepPhi

    SUBROUTINE GaussSeidel(data_in)
        !> Subroutine to perform the Gauss-Seidel algorithm to a Poisson equation.
        TYPE(GaussSeidelData), INTENT(INOUT)    :: data_in
        INTEGER                                 :: t_idx
        LOGICAL                                 :: is_converged

        is_converged = .FALSE.

        DO t_idx = 1, data_in%n_iter

            DO WHILE ( .NOT. is_converged)

                CALL SweepPhi(data_in)
                is_converged = ConvergenceTest(data_in)

            END DO

            !> Update E fields once suitable phi convergence is attained.
            CALL EfieldUpdate(data_in)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Insert velocity-verlet stuff here !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        END DO

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Insert netCDF stuff here !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END SUBROUTINE GaussSeidel

    SUBROUTINE Display_Help()
        !> Provides a help option for users unsure of correct inputs.
        !! Command line arguments are stated with default values listed if argument is missing.

        PRINT *, "                                                                                                  "
        PRINT *, "                        Electrostatics [options]                                                  "
        PRINT *, "  Options:                                                                       Default Setting: "
        PRINT *, "                                                                                                  "
        PRINT *, "  problem = <char>        Initialisation type: null, single, double                    -          "
        PRINT *, "                          -null: rho = 0                                                          "
        PRINT *, "                          -single: rho = EXP(-(x/0.1)^2 - (y/0.1)^2)                              "
        PRINT *, "                          -double: rho = EXP(-((x+0.25)/0.1)^2 - ((y+0.25)/0.1)^2)                &
              &                                          + EXP(-((x-0.75)/0.2)^2 - ((y-0.75)/0.2)^2)                "
        PRINT *, "  nx = <integer>          X-grid length                                                50         "
        PRINT *, "  ny = <integer>          Y-grid length                                                50         "
        PRINT *, "  n_iter = <integer>      Number of iterations                                        1000        "
        PRINT *, "  dt = <real64>           Time steps                                                  0.01        "

    RETURN
    END SUBROUTINE Display_Help

    SUBROUTINE Parse_Command_Line(init_type, nx, ny, n_iter, dt)

        CHARACTER(LEN=*), INTENT(OUT)       :: init_type
        INTEGER, INTENT(OUT)                :: nx, ny
        INTEGER, INTENT(OUT), OPTIONAL      :: n_iter
        REAl(REAL64), INTENT(OUT), OPTIONAL :: dt

        INTEGER                             :: arg_count, arg_idx, temp_int
        REAL(REAL64)                        :: temp_real
        LOGICAL                             :: Exists, IsValid
        CHARACTER(LEN=20)                   :: temp_str, arg

        arg_count = COMMAND_ARGUMENT_COUNT()

        DO arg_idx = 1, arg_count
            CALL GET_COMMAND_ARGUMENT(arg_idx, arg)
            arg = TRIM(arg)
            IF (arg == '--help') THEN
                !> Custom help command for command line inputs 
                CALL Display_Help()
                STOP

            END IF
        END DO

        CALL parse_args()

        IsValid = get_arg('problem', temp_str, Exists)
        IF (IsValid .AND. Exists) THEN
            SELECT CASE(temp_str)
            CASE('null')
                init_type = 'null'

            CASE('single')
                init_type = 'single'

            CASE('double')
                init_type = 'double'

            CASE DEFAULT
                CALL Add_Error_Message("Error: Invalid initialisation selected.")
                CALL Print_Errors()
                ERROR STOP

            END SELECT


        END IF


        IsValid = get_arg('nx', temp_int, Exists)
        IF (IsValid .AND. Exists) THEN
            IF (temp_int > 0) THEN
                !> Assigns grid length to correct variable.
                nx = temp_int

            ELSE
                !> Warn user that default variable is being used due to invalid grid size.
                !! Adds corresponding error message to error log.
                CALL Add_Error_Message("Error: Grid length nx must be a positive integer.")
                CALL Print_Errors()
                ERROR STOP

            END IF

        ELSE
            !> Terminate program due to missing argument.
            CALL Add_Error_Message("Error: Grid length nx not provided.")
            CALL Print_Errors()
            ERROR STOP
        END IF

        IsValid = get_arg('ny', temp_int, Exists)
        IF (IsValid .AND. Exists) THEN
            IF (temp_int > 0) THEN
                !> Assigns grid length to correct variable.
                ny = temp_int

            ELSE
                !> Warn user that default variable is being used due to invalid grid size.
                !! Adds corresponding error message to error log.
                CALL Add_Error_Message("Error: Grid length ny must be a positive integer.")
                CALL Print_Errors()
                ERROR STOP

            END IF

        ELSE
            !> Terminate program due to missing argument.
            CALL Add_Error_Message("Error: Grid length ny not provided.")
            CALL Print_Errors()
            ERROR STOP
        END IF

        IsValid = get_arg('n_iter', temp_int, Exists)
        IF (IsValid .AND. Exists) THEN
            IF (temp_int > 0) THEN
                !> Assigns grid length to correct variable.
                n_iter = temp_int

            ELSE
                !> Warn user that default variable is being used due to invalid grid size.
                !! Adds corresponding error message to error log.
                CALL Add_Error_Message("Error: Number of iterations must be a positive integer.")
                CALL Add_Warning_Message("Warning: Using default value for number of iterations .")
                n_iter = 1000

            END IF

        ELSE
            !> Set to default setting.
            CALL Add_Warning_Message("Warning: Number of iterations not provided. Using default value.")
            n_iter = 1000
        END IF

        IsValid = get_arg('dt', temp_real, Exists)
        IF (IsValid .AND. Exists) THEN
            IF (temp_int > 0) THEN
                !> Assigns grid length to correct variable.
                dt = temp_real

            ELSE
                !> Warn user that default variable is being used due to invalid grid size.
                !! Adds corresponding error message to error log.
                CALL Add_Error_Message("Error: Time step must be a positive real.")
                CALL Add_Warning_Message("Warning: Using default value for time step.")
                dt = 0.01_REAL64

            END IF

        ELSE
            !> Terminate program due to missing argument.
            CALL Add_Warning_Message("Warning: Time step not provided. Using default value.")
            dt = 0.01_REAL64
        END IF

    END SUBROUTINE Parse_Command_Line

END MODULE solve_gauss_seidel
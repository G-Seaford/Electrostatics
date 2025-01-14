MODULE SolveGaussSeidel
    !> Solves the Poisson equation iteratively via the Gauss-Seidel method.
    !! Uses the Velocity Verlet algorithm by Facundo Costa to update positions, accelerations and velocities.
    !! Writes data to a NetCDF file for plotting in external software.
    !! @author Gianluca Seaford
    !! @version 2.0
    !! @date 2024-11-29

    USE command_line
    USE CommonDataStructure
    USE Error_Logging
    USE ISO_FORTRAN_ENV
    USE VelocityVerlet
    USE write_netcdf

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE InitialiseCommonData(data, status)
        !> Initialises the data structure for the Gauss-Seidel method.
        !! @param[inout] data CommonData structure to initialise.
        !! @param[inout] status Status for error handling.

        TYPE(CommonData), INTENT(INOUT)     :: data
        INTEGER, INTENT(INOUT)              :: status

        CHARACTER(LEN=6)                    :: init_type
        INTEGER                             :: nx, ny, n_iter, i, j, cell_x, cell_y
        REAL(REAL64)                        :: x, y, dt
        

        CALL ParseCommandLine(init_type, nx, ny, n_iter, dt)

        !> Set parameters within CommonData structure.
        data%nx = nx
        data%ny = ny
        data%n_iter = n_iter

        data%dx = 2.0 / nx
        data%dy = 2.0 / ny
        data%dt = dt

        !> Allocate field arrays.
        !! If an error is detected, the CleanUpData subroutine
        !! is called to deallocate allocated arrays.

        ALLOCATE(data%phi(0:nx+1, 0:ny+1), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for phi")
            CALL CleanUpData(data)
            ERROR STOP
        END IF

        ALLOCATE(data%rho(1:nx, 1:ny), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for rho")
            CALL CleanUpData(data)
            ERROR STOP
        END IF

        ALLOCATE(data%Ex(1:nx, 1:ny), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for Ex")
            CALL CleanUpData(data)
            ERROR STOP
        END IF

        ALLOCATE(data%Ey(1:nx, 1:ny), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for Ey")
            CALL CleanUpData(data)
            ERROR STOP
        END IF

        !> Allocate particle storage arrays.
        !! If an error is detected, the CleanUpData subroutine
        !! is called to deallocate allocated arrays.

        ALLOCATE(data%positions(1:2, 0:n_iter), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for positions")
            CALL CleanUpData(data)
            ERROR STOP
        END IF

        ALLOCATE(data%velocities(1:2, 0:n_iter), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for velocities")
            CALL CleanUpData(data)
            ERROR STOP
        END IF
        
        ALLOCATE(data%accelerations(1:2, 0:n_iter), STAT=status)
        IF (status /= 0) THEN
            CALL Add_Error_Message("Error: Failed allocating memory for accelerations")
            CALL CleanUpData(data)
            ERROR STOP
        END IF

        !> Initialise arrays with zero values.
        !! This is equivalent to imposing Dirichlet conditions at boundaries.
        data%phi = 0.0_REAL64
        data%rho = 0.0_REAL64

        data%Ex = 0.0_REAL64
        data%Ey = 0.0_REAL64

        data%positions = 0.0_REAL64
        data%velocities = 0.0_REAL64
        data%accelerations = 0.0_REAL64

        SELECT CASE(init_type)

            CASE("null")
                !> Zero Gaussian peaks, with null rho.
                !! Used for a lone charge.

                !> Set initial x- and y-velocity.
                data%velocities(1:2, 0) = 0.1_REAL64

                !> Set initial E-fields.
                CALL EfieldUpdate(data)

                !> Set initial acceleration
                cell_x = FLOOR((data%positions(1,0) + 1.0_REAL64)/data%dx) + 1
                cell_y = FLOOR((data%positions(2,0) + 1.0_REAL64)/data%dy) + 1

                data%accelerations(1,0) = -1.0 * data%Ex(cell_x,cell_y)
                data%accelerations(2,0) = -1.0 * data%Ey(cell_x,cell_y)
                
                RETURN

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
                
                !> Set initial E-fields.
                CALL EfieldUpdate(data)

                !> Set initial x velocity.
                data%positions(1,0) = 0.1_REAL64

                !> Set initial acceleration
                cell_x = FLOOR((data%positions(1,0) + 1.0_REAL64)/data%dx) + 1
                cell_y = FLOOR((data%positions(2,0) + 1.0_REAL64)/data%dy) + 1

                data%accelerations(1,0) = -1.0 * data%Ex(cell_x,cell_y)
                data%accelerations(2,0) = -1.0 * data%Ey(cell_x,cell_y)

                RETURN

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

                !> Set initial E-fields.
                CALL EfieldUpdate(data)

                !> Set initial y position.
                data%positions(2,0) = 0.5_REAL64

                !> Set initial acceleration
                cell_x = FLOOR((data%positions(1,0) + 1.0_REAL64)/data%dx) + 1
                cell_y = FLOOR((data%positions(2,0) + 1.0_REAL64)/data%dy) + 1

                data%accelerations(1,0) = -1.0 * data%Ex(cell_x,cell_y)
                data%accelerations(2,0) = -1.0 * data%Ey(cell_x,cell_y)

                RETURN

            CASE DEFAULT
                CALL Add_Error_Message("Error: Unknown initialisation type.")
                CALL Print_Errors()

                ERROR STOP

        END SELECT

    END SUBROUTINE InitialiseCommonData

    SUBROUTINE CleanUpData(data)
        TYPE(CommonData), INTENT(INOUT) :: data

        INTEGER                         :: status

        status = 0

        !> Deallocate field arrays.
        IF (ALLOCATED(data%phi)) THEN
            DEALLOCATE(data%phi, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for phi")
                RETURN
            END IF
        END IF
        
        IF (ALLOCATED(data%rho)) THEN
            DEALLOCATE(data%rho, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for rho")
                RETURN
            END IF
        END IF

        IF (ALLOCATED(data%Ex)) THEN
            DEALLOCATE(data%Ex, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for Ex")
                RETURN
            END IF
        END IF

        IF (ALLOCATED(data%Ey)) THEN
            DEALLOCATE(data%Ey, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for Ey")
                RETURN
            END IF
        END IF

        !> Deallocate particle storage arrays.
        IF (ALLOCATED(data%positions)) THEN
            DEALLOCATE(data%positions, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for positions")
                RETURN
            END IF
        END IF

        IF (ALLOCATED(data%velocities)) THEN
            DEALLOCATE(data%velocities, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for velocities")
                RETURN
            END IF
        END IF
        
        IF (ALLOCATED(data%accelerations)) THEN
            DEALLOCATE(data%accelerations, STAT=status)
            IF (status /= 0) THEN
                CALL Add_Error_Message("Error: Failed deallocating memory for accelerations")
                RETURN
            END IF
        END IF

        CALL Print_Errors()

    END SUBROUTINE CleanUpData

    SUBROUTINE PhiUpdate(data_in, x_idx, y_idx)
        !> Subroutine to update a given phi element.
        !! Calculates the second-order partial derivatives via a central finite difference method.
        !! @param[inout] data_in CommonData structure holding data.
        !! @param[in] x_idx The x index.
        !! @param[in] y_idx The y index.

        TYPE(CommonData), INTENT(INOUT)     :: data_in
        INTEGER, INTENT(IN)                 :: x_idx, y_idx

        REAL(REAL64)                        :: dx, dy, dx2_inv, dy2_inv, rho, phi_current, phi_xdiff, phi_ydiff, denom_inv

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
        !! @param[inout] data_in CommonData structure holding data.

        TYPE(CommonData), INTENT(INOUT)     :: data_in

        INTEGER                             :: x_idx, y_idx
        REAL(REAL64)                        :: dx, dy, phi_xdiff, phi_ydiff

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
        !! @param[in] data_in CommonData structure holding data.
        !! @param[out] isconverged Logical determining whether convergence has been attained.

        TYPE(CommonData), INTENT(IN)    :: data_in

        INTEGER                         :: x_idx, y_idx, n_sites
        REAL(REAL64)                    :: rho, dx, dy, dx2_inv, dy2_inv, phi_xdiff, phi_ydiff, tot_diff, &
                                        &  e_tot, d_rms, ratio, epsilon = 1.0e-12_REAL64

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
            !> Prevent singularities by not dividing by suitably small d_rms values.

            IF (e_tot < epsilon) THEN
                !> True convergence achieved.
                !! Because both the total error and RMS error ers suitably small.
                ratio = 0.0_REAL64
            
            ELSE
                !> Indicative of an ill-conditioned Laplacian matrix.
                !! Since the total error is significantly larger than the RMS error.
                ratio = HUGE(1.0_REAL64)
            END IF 

        ELSE 
            !> Evaluate the ratio of the total error to the RMS distance
            ratio = e_tot / d_rms

        END IF

        IF (ratio < 0.00001_REAL64) THEN
            !> Check convergence criterion.
            !! e_tot/d_rms < 10^-5.
            is_converged = .TRUE.
        END IF

    END FUNCTION ConvergenceTest

    SUBROUTINE SweepPhi(data_in)
        !> Performs a Gauss–Seidel sweep over the entire phi array,
        !! updating phi at each grid point.
        !! @param[inout] data_in CommonData structure containing the grid and field arrays.

        TYPE(CommonData), INTENT(INOUT)     :: data_in
        INTEGER                             :: x_idx, y_idx

        DO y_idx=1, data_in%ny
            DO x_idx=1, data_in%nx

                CALL PhiUpdate(data_in, x_idx, y_idx)

            END DO

        END DO

    END SUBROUTINE SweepPhi

    SUBROUTINE VelocityVerletWrapper(data_in)
        !> Wrapper to call the velocity verlet algorithm.
        !! @param[inout] data_in CommonData structure.

        TYPE(CommonData), INTENT(INOUT)     :: data_in

        CALL velocity_verlet(data_in%positions, data_in%velocities, data_in%accelerations, & 
                           & data_in%Ex, data_in%Ey, data_in%n_iter, data_in %nx, data_in%ny, &
                           & data_in%dx, data_in%dy, data_in%dt)

    END SUBROUTINE VelocityVerletWrapper

    SUBROUTINE NetCDFWriteWrapper(data_in, fname_in)
        !> Wrapper to call the NetCDF writing function.
        !! @param[in] data_in CommonData structure.
        !! @param[in] fname_in Name for the netCDF file.

        TYPE(CommonData), INTENT(IN)    :: data_in
        CHARACTER(len=*), INTENT(IN)    :: fname_in

        INTEGER                         :: idx_err

        CALL writer(data_in, fname_in, idx_err)

    END SUBROUTINE NetCDFWriteWrapper

    SUBROUTINE SolveSystem(data_in)
        !> Subroutine to perform the Gauss-Seidel algorithm to a Poisson equation.
        TYPE(CommonData), INTENT(INOUT)     :: data_in

        INTEGER                             :: sweep_idx
        LOGICAL                             :: is_converged

        is_converged = .FALSE.
        sweep_idx = 0

        DO WHILE ( .NOT. is_converged .AND. sweep_idx < 1000)

            CALL SweepPhi(data_in)
            is_converged = ConvergenceTest(data_in)
            sweep_idx = sweep_idx + 1

        END DO

        !> Update E fields once suitable phi convergence is attained.
        CALL EfieldUpdate(data_in)

        CALL VelocityVerletWrapper(data_in)

        CALL NetCDFWriteWrapper(data_in, 'Gauss-Seidel_Electrostatics.nc')

    END SUBROUTINE SolveSystem

    SUBROUTINE DisplayHelp()
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
        PRINT *, "                                                                                                  "

    RETURN

    END SUBROUTINE DisplayHelp

    SUBROUTINE ParseCommandLine(init_type, nx, ny, n_iter, dt)

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
                CALL DisplayHelp()
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
            CALL DisplayHelp()
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
                CALL DisplayHelp()
                CALL Print_Errors()
                ERROR STOP

            END IF

        ELSE
            !> Terminate program due to missing argument.
            CALL Add_Error_Message("Error: Grid length ny not provided.")
            CALL DisplayHelp()
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
            IF (temp_real > 0.0_REAL64) THEN
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

    END SUBROUTINE ParseCommandLine

END MODULE SolveGaussSeidel
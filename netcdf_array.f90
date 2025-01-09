MODULE write_netcdf

  USE ISO_FORTRAN_ENV
  USE netcdf
  USE run_data_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE writer(run, filename, ierr)
    TYPE(run_data), INTENT(IN) :: run   ! Pass the run_data type containing all necessary variables
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(INOUT) :: ierr
    INTEGER, PARAMETER :: ndims2 = 2, ndims1 = 1
    INTEGER, DIMENSION(ndims2) :: sizes2, dim_ids2
    INTEGER, DIMENSION(ndims1) :: sizes1, dim_ids1
    INTEGER :: file_id, var_id_2d, var_id_3d, var_id_1d, i
    INTEGER :: size_t, size_ex, size_ey, size_potential
    CHARACTER(LEN=1), DIMENSION(ndims2) :: dims2 = (/"x", "y"/)
    CHARACTER(LEN=1), DIMENSION(ndims1) :: dims1 = (/"t"/)

    ! Sizes of the arrays
    sizes2 = SHAPE(run%position)
    sizes1 = SHAPE(run%axis_x)
    size_t = SIZE(run%axis_x)

    ! Get the dimensions of the potential, E_x, and E_y
    size_ex = SIZE(run%E_x)
    size_ey = SIZE(run%E_y)
    size_potential = SIZE(run%potential)
    
    ! Debug: Print the sizes of arrays before writing
    PRINT*, "Sizes of position (2D array): ", sizes2
    PRINT*, "Size of axis_x (1D array): ", sizes1
    PRINT*, "Size of E_x (3D array): ", size_ex
    PRINT*, "Size of E_y (3D array): ", size_ey
    PRINT*, "Size of potential (3D array): ", size_potential

    ! Create the file, overwriting if it exists, because this can be run multiple times
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error creating NetCDF file: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "File created successfully."

    ! Define dimensions for the rank-2 (2D arrays like position, velocity, etc.)
    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i), sizes2(i), dim_ids2(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    PRINT*, "Dimensions for position, velocity, acceleration, etc. defined successfully."

    ! Define dimensions for the rank-1 (axis_x and axis_y)
    DO i = 1, ndims1
      ierr = nf90_def_dim(file_id, dims1(i), sizes1(i), dim_ids1(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for ", dims1(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    PRINT*, "Rank-1 dimensions (axis_x, axis_y) defined successfully."

    ! Define the rank-2 variables (e.g., position, velocity, etc.)
    ierr = nf90_def_var(file_id, "position", NF90_REAL, dim_ids2, var_id_2d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'position': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'position' variable defined successfully."

    ierr = nf90_def_var(file_id, "velocity", NF90_REAL, dim_ids2, var_id_3d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'velocity': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'velocity' variable defined successfully."

    ierr = nf90_def_var(file_id, "acceleration", NF90_REAL, dim_ids2, var_id_3d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'acceleration': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'acceleration' variable defined successfully."

    ! Define dimensions for potential, E_x, E_y based on their actual shapes
    ! For example, if these are 3D arrays, we'll need to define dimensions for each
    ! Assuming their shapes are of the form (nx, ny, nt) or something similar

    ! Define dimensions for potential, E_x, and E_y (dynamically based on their shape)
    ! Here we're using size_t, size_ex, size_ey, and size_potential to handle the shape
    ierr = nf90_def_dim(file_id, "nx", size_ex, dim_ids2(1))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining dimension for 'nx': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_dim(file_id, "ny", size_ey, dim_ids2(2))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining dimension for 'ny': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Define potential, E_x, and E_y variables using these dimensions
    ierr = nf90_def_var(file_id, "potential", NF90_REAL, dim_ids2, var_id_3d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'potential': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'potential' variable defined successfully."

    ierr = nf90_def_var(file_id, "E_x", NF90_REAL, dim_ids2, var_id_3d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'E_x': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_x' variable defined successfully."

    ierr = nf90_def_var(file_id, "E_y", NF90_REAL, dim_ids2, var_id_3d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'E_y': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_y' variable defined successfully."

    ! Define the rank-1 variables for axis_x and axis_y
    ierr = nf90_def_var(file_id, "axis_x", NF90_REAL, dim_ids1, var_id_1d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'axis_x': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'axis_x' variable defined successfully."

    ierr = nf90_def_var(file_id, "axis_y", NF90_REAL, dim_ids1, var_id_1d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'axis_y': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'axis_y' variable defined successfully."

    ! Define attributes for the global data (e.g., num_steps, nx, ny)
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "num_steps", run%num_steps)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error setting global attribute 'num_steps': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "Global attribute 'num_steps' set successfully."

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "nx", run%nx)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error setting global attribute 'nx': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "Global attribute 'nx' set successfully."

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "ny", run%ny)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error setting global attribute 'ny': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "Global attribute 'ny' set successfully."

    ! Finish defining metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error finishing NetCDF definitions: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "NetCDF definitions finished successfully."

    ! Write the rank-2 variables (e.g., position, velocity, etc.)
    ierr = nf90_put_var(file_id, var_id_2d, run%position)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'position' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'position' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_3d, run%velocity)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'velocity' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'velocity' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_3d, run%acceleration)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'acceleration' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'acceleration' variable written successfully."

    ! Write the potential, E_x, and E_y variables
    ierr = nf90_put_var(file_id, var_id_3d, run%potential)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'potential' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'potential' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_3d, run%E_x)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'E_x' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_x' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_3d, run%E_y)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'E_y' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_y' variable written successfully."

    ! Write the rank-1 variables (axis_x and axis_y)
    ierr = nf90_put_var(file_id, var_id_1d, run%axis_x)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'axis_x' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'axis_x' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_1d, run%axis_y)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'axis_y' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'axis_y' variable written successfully."

    ! Close
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error closing NetCDF file: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "NetCDF file closed successfully."

  END SUBROUTINE writer

END MODULE write_netcdf

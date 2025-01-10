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
    INTEGER, PARAMETER :: ndims2 = 2
    INTEGER, DIMENSION(ndims2) :: sizes2, dim_ids2, sizes_ex, sizes_ey, sizes_potential, dim_ids_ex, dim_ids_ey, dim_ids_potential
    INTEGER :: file_id, var_id_2d, var_id_ex, var_id_ey, var_id_potential, i
    CHARACTER(LEN=1), DIMENSION(ndims2) :: dims2 = (/'x', 'y'/)

    ! Sizes of the arrays
    sizes2 = SHAPE(run%position)
    sizes_ex = SHAPE(run%E_x)
    sizes_ey = SHAPE(run%E_y)
    sizes_potential = SHAPE(run%potential)

    ! Debug: Print the sizes of arrays before writing
    PRINT*, "Sizes of position (2D array): ", sizes2
    PRINT*, "Sizes of E_x (2D array): ", sizes_ex
    PRINT*, "Sizes of E_y (2D array): ", sizes_ey
    PRINT*, "Sizes of potential (2D array): ", sizes_potential

    ! Create the file, overwriting if it exists, because this can be run multiple times
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error creating NetCDF file: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "File created successfully."

    ! Define dimensions for the rank-2 (2D arrays like position, E_x, E_y, potential)
    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i), sizes2(i), dim_ids2(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    PRINT*, "Dimensions for position defined successfully."

    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i) // "_ex", sizes_ex(i), dim_ids_ex(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for E_x ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    PRINT*, "Dimensions for E_x defined successfully."

    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i) // "_ey", sizes_ey(i), dim_ids_ey(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for E_y ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    PRINT*, "Dimensions for E_y defined successfully."

    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i) // "_pot", sizes_potential(i), dim_ids_potential(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for potential ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    PRINT*, "Dimensions for potential defined successfully."

    ! Define the rank-2 variables
    ierr = nf90_def_var(file_id, "position", NF90_REAL, dim_ids2, var_id_2d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'position': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'position' variable defined successfully."
    
    ierr = nf90_def_var(file_id, "velocity", NF90_REAL, dim_ids2, var_id_2d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable velocity'': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'velocity' variable defined successfully."

    ierr = nf90_def_var(file_id, "acceleration", NF90_REAL, dim_ids2, var_id_2d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable acceleration'': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'acceleration' variable defined successfully."


    ierr = nf90_def_var(file_id, "E_x", NF90_REAL, dim_ids_ex, var_id_ex)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'E_x': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_x' variable defined successfully."

    ierr = nf90_def_var(file_id, "rho", NF90_REAL, dim_ids_ex, var_id_ex)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'rho': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'rho' variable defined successfully."

    ierr = nf90_def_var(file_id, "E_y", NF90_REAL, dim_ids_ey, var_id_ey)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'E_y': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_y' variable defined successfully."

    ierr = nf90_def_var(file_id, "potential", NF90_REAL, dim_ids_potential, var_id_potential)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'potential': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'potential' variable defined successfully."

    ! Finish defining metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error finishing NetCDF definitions: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "NetCDF definitions finished successfully."

    ! Write the rank-2 variables
    ierr = nf90_put_var(file_id, var_id_2d, run%position)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'position' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'position' variable written successfully."

    ! Write the rank-2 variables
    ierr = nf90_put_var(file_id, var_id_2d, run%acceleration)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'acceleration' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'acceleration' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_ex, run%E_x)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'E_x' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_x' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_ex, run%rho)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'rho' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'rho' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_ey, run%E_y)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'E_y' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'E_y' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_potential, run%potential)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'potential' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'potential' variable written successfully."

    ! Close the NetCDF file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error closing NetCDF file: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "NetCDF file closed successfully."

  END SUBROUTINE writer

END MODULE write_netcdf

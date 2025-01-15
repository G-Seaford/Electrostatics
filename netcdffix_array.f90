MODULE write_netcdf
  ! Module to write data to a netCDF file, which when will then read after
  ! by facu
  USE CommonDataStructure
  USE ISO_FORTRAN_ENV
  USE netcdf

  IMPLICIT NONE

CONTAINS

  SUBROUTINE writer(run, filename, ierr)
    TYPE(CommonData), INTENT(IN) :: run   ! Pass the run_data type containing all necessary variables
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(INOUT) :: ierr
    INTEGER, PARAMETER :: ndims2 = 2
    INTEGER, DIMENSION(ndims2) :: sizes2, dim_ids2, sizes_ex, sizes_ey, sizes_phi, dim_ids_ex, dim_ids_ey, dim_ids_phi
    INTEGER :: file_id, var_id_2d, var_id_ex, var_id_ey, var_id_phi, i,var_id_vel,var_id_acc,var_id_rho
    CHARACTER(LEN=1), DIMENSION(ndims2) :: dims2 = (/'x', 'y'/)

    ! Sizes of the arrays
    sizes2 = SHAPE(run%positions)
    sizes_ex = SHAPE(run%Ex)
    sizes_ey = SHAPE(run%Ey)
    sizes_phi = SHAPE(run%phi)

    ! Debugging so Print the sizes of arrays before writing
    !PRINT*, "Sizes of position (2D array): ", sizes2
    !PRINT*, "Sizes of E_x (2D array): ", sizes_ex
    !PRINT*, "Sizes of E_y (2D array): ", sizes_ey
    !PRINT*, "Sizes of phi (2D array): ", sizes_phi

    ! Create the file, overwriting if it exists, because this can be run multiple times
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error creating NetCDF file: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "File created successfully."

    ! Define dimensions for accel,position,velo
    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i), sizes2(i), dim_ids2(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i) // "_ex", sizes_ex(i), dim_ids_ex(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for Ex ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO
    

    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i) // "_ey", sizes_ey(i), dim_ids_ey(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for Ey ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    DO i = 1, ndims2
      ierr = nf90_def_dim(file_id, dims2(i) // "_pot", sizes_phi(i), dim_ids_phi(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, "Error defining dimension for phi ", dims2(i), ": ", TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ! I stopped defining the dimension for every variable as lots have the sames ones

    ! Define the rank-2 variables

    ! rho has same dim as E_x so just re-used this, else I would not do this
    ierr = nf90_def_var(file_id, "rho", NF90_REAL, dim_ids_ex, var_id_rho)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'rho': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'rho' variable defined successfully."

    
    ierr = nf90_def_var(file_id, "phi", NF90_REAL, dim_ids_phi, var_id_phi)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'phi': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'phi' variable defined successfully."

    ierr = nf90_def_var(file_id, "Ex", NF90_REAL, dim_ids_ex, var_id_ex)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'Ex': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'Ex' variable defined successfully."

    ierr = nf90_def_var(file_id, "Ey", NF90_REAL, dim_ids_ey, var_id_ey)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'Ey': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'Ey' variable defined successfully."

    ierr = nf90_def_var(file_id, "positions", NF90_REAL, dim_ids2, var_id_2d)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable 'positions': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'positions' variable defined successfully."
    
    ierr = nf90_def_var(file_id, "velocities", NF90_REAL, dim_ids2, var_id_vel)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable velocities'': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'velocities' variable defined successfully."

    ierr = nf90_def_var(file_id, "accelerations", NF90_REAL, dim_ids2, var_id_acc)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error defining variable accelerations'': ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'accelerations' variable defined successfully."

    ! Finish defining metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error finishing NetCDF definitions: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "NetCDF definitions finished successfully."

    ! Write the rank-2 variables
    ierr = nf90_put_var(file_id, var_id_rho, run%rho)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'rho' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'rho' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_phi, run%phi)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'phi' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'phi' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_ex, run%Ex)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'Ex' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'Ex' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_ey, run%Ey)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'Ey' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'Ey' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_2d, run%positions)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'positions' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'positions' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_vel, run%velocities)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'velocities' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'velocities' variable written successfully."

    ierr = nf90_put_var(file_id, var_id_acc, run%accelerations)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error writing 'accelerations' variable: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "'accelerations' variable written successfully."

    ! Close the NetCDF file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, "Error closing NetCDF file: ", TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    PRINT*, "NetCDF file closed successfully."

  END SUBROUTINE writer

END MODULE write_netcdf

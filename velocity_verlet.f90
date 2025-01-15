MODULE VelocityVerlet
  !> A module that moves the electron which is a test charge around, which uses the velocity verlet, and the force which is coming from grad of the potential at the points.
  ! by facu

  ! first, find out where on earth the electron is in terms of the grid
  ! second thing, from the scalar potentail, take the finite difference in the x and the y and then from this, get the electric field at the point in the grid.
  
  ! We used lucas E_field solver instead but they both agreed
  
  ! using the velocity verlet algo, find the new position based on the old one, the old velocity and old accel
  ! then find the new accelerations
  ! and the new velocity
  ! and then loop round to finding the new position and so on....
  ! after the field is set up, we need to run the velocity verlet for 1000 time steps that is the idea
  ! vars: current_x, current_y, current_vx, current_vy, current_ax, current_ay
  
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  SAVE
  CONTAINS

  SUBROUTINE velocity_verlet(position, velocity, acceleration, E_field_x, E_field_y, num_steps, nx, ny, dx, dy, delta_t)
    ! charge and mass are set to 1 and so is the permitivity of free space
    ! Input/Output Variables

    INTEGER, INTENT(IN) :: num_steps, nx, ny
    REAL(REAL64), INTENT(IN) :: delta_t, dx, dy
    REAL(REAL64), DIMENSION(:,:), INTENT(IN) :: E_field_x, E_field_y
    REAL(REAL64), DIMENSION(1:2,0:num_steps) :: position, velocity, acceleration

    INTEGER :: cell_x, cell_y, step, endsteps

    ! Main Velocity Verlet Loop

    DO step = 0, num_steps - 1
  
      ! updating positions

      position(1,step+1) = position(1,step) + velocity(1,step) * delta_t + 0.5_REAL64 * acceleration(1,step) * delta_t**2
      position(2,step+1) = position(2,step) + velocity(2,step) * delta_t + 0.5_REAL64 * acceleration(2,step) * delta_t**2 
      

      ! finds what cell the new position is attributed too
      cell_x = FLOOR((position(1,step+1) + 1.0_REAL64)/dx) + 1
      cell_y = FLOOR((position(2,step+1) + 1.0_REAL64)/dy) + 1

      
      ! the below deals with the electron having left the grid and stops the program
      ! where all later positions/velocties/accelerations are the same as the one when
      ! it left, as we stop the algorithm

      if ((cell_x < 1 .or. cell_x > nx) .or. (cell_y < 1 .or. cell_y > ny)) then


        ! This will fill all the data left
        DO endsteps = step + 1, num_steps
          position(1,endsteps) = position(1,step)
          position(2,endsteps) = position(2,step)

          acceleration(1,endsteps) = acceleration(1,step)
          acceleration(2,endsteps) = acceleration(2,step)

          velocity(1,endsteps) =  velocity(1,step)
          velocity(2,endsteps) =   velocity(2,step)
        end DO
      else

        ! find acceleration using newtons 2nd law

        acceleration(1,step+1) = -1.0*E_field_x(cell_x,cell_y)
        acceleration(2,step+1) = -1.0*E_field_y(cell_x,cell_y)

        ! update the velocity
        velocity(1,step+1) = velocity(1,step) + 0.5_REAL64*delta_t*(acceleration(1,step + 1)+acceleration(1,step))
        velocity(2,step+1) = velocity(2,step) + 0.5_REAL64*delta_t*(acceleration(2,step + 1)+acceleration(2,step))
      end if

    END DO

END SUBROUTINE velocity_verlet

END MODULE VelocityVerlet

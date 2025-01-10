PROGRAM main
  USE domain_tools
  USE move_elec
  USE ISO_FORTRAN_ENV
  USE run_data_module
  USE write_netcdf
  IMPLICIT NONE

  TYPE(run_data) :: run
  REAL(REAL64), DIMENSION(:), ALLOCATABLE :: axis_x
  REAL(REAL64), DIMENSION(:), ALLOCATABLE :: axis_y
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: E_x
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: E_y
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: position
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: velocity
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: acceleration
  INTEGER :: num_steps

  INTEGER :: nx
  INTEGER :: ny
  REAL(REAL64), DIMENSION(2) :: axis_range_x
  REAL(REAL64), DIMENSION(2) :: axis_range_y
  INTEGER :: nghosts
  INTEGER :: err
  REAL(REAL64) ::dx
  REAL(REAL64) ::dy
  REAL(REAL64) ::delta_t
  REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: potential

  num_steps = 1000
  delta_t = 0.01
  nx = 100
  ny = 100
  axis_range_x = [-1.0_REAL64, 1.0_REAL64]
  axis_range_y = [-1.0_REAL64, 1.0_REAL64]
  nghosts = 1

  CALL create_axis(axis_x, nx, axis_range_x, nghosts,dx)
  CALL create_axis(axis_y, ny, axis_range_y, nghosts,dy)

  ! print*, axis_x
  ! print*, axis_y

  ! grid is all working !


  ! here we have the grid created

  ALLOCATE(potential(0:nx+1, 0:ny+1))
  ALLOCATE(position(2,num_steps))
  ALLOCATE(velocity(2,num_steps))
  ALLOCATE(acceleration(2,num_steps))




  position(1,1) = 0.0_REAL64
  position(2,1) = 0.0_REAL64
  velocity(1,1) = -0.01_REAL64
  velocity(2,1) = 0.002_REAL64
  !Then run gauss siedel here to get what the potentail should be!

  potential = 1.0_REAL64



  !Then run create E_field from the potentail




  CALL create_E_field(potential, nx, ny, dx, dy, E_x, E_y)
  
  ! Now that we have the E field, we can run the simulation
  CALL velocity_verlet(position, velocity, acceleration, E_x, E_y, delta_t,num_steps,dx,dy,nx,ny)

  
  ! USE NETcdf now so that we can ouput all the figures nicely

  run%axis_x = axis_x
  run%axis_y = axis_y
  run%E_x = E_x
  run%E_y = E_y
  run%position = position
  run%velocity = velocity
  run%acceleration = acceleration
  run%potential = potential
  run%num_steps = num_steps
  run%nx = nx
  run%ny = ny
  

  CALL writer(run,"saved.nc", err)
  if(err /= 0)then
    print*, "Error writing to the file"
  end if


  DEALLOCATE(position)
  DEALLOCATE(velocity)
  DEALLOCATE(acceleration)
  DEALLOCATE(potential)
  DEALLOCATE(E_x)
  DEALLOCATE(E_y)
  DEALLOCATE(axis_x)
  DEALLOCATE(axis_y)

END PROGRAM main

program surfacehopping
  use dynamik
  implicit none
  
  integer :: trajectory
  ! parameters read from input file
  integer :: N_traj
  integer :: state_init
  double precision :: r_init, p_init
  double precision :: mass
  double precision :: dE, angle_dp
  double precision :: upper_limitZ, lower_limitZ
  double precision :: dt
  ! variables for statistics
  integer :: regular_traj, irregular_traj, end_state(N_states)
  double precision :: deltaE_max, deltaNorm_max, population(N_states)


  open(100,file='trajinfo.out',status='replace')

  ! read input parameters
  call read_input_parameters(N_traj, state_init, r_init, p_init, mass, dE, upper_limitZ, lower_limitZ, dt)

  ! read potentials and couplings
  call read_input

  ! set random seed
  call set_nran
  
  ! initialize statistics variables
  regular_traj = 0
  end_state(:) = 0
  population(:) = 0.d0
  deltaE_max = 0.d0
  deltaNorm_max = 0.d0
  call init_histo()

  all_trajs: do trajectory = 1,N_traj
     
     call set_initial_cond(state_init, r_init, p_init, mass, dE, dt)
     single_traj: do

        call timestep
        !projectile is reflected from the surface
        if (ret_actual_posZ().gt.upper_limitZ) then
           regular_traj = regular_traj + 1
           end_state(ret_actual_state()) = end_state(ret_actual_state()) + 1
           population(:) = population(:) + ret_actual_population()
           call fill_histo()
           exit
        end if
        !projectile enters the solid
        if (ret_actual_posZ().lt.lower_limitZ) exit

     end do single_traj

     call print_traj_info()
     if (deltaE_max.lt.ret_deltaE_max()) deltaE_max = ret_deltaE_max()
     if (deltaNorm_max.lt.ret_deltaNorm_max()) deltaNorm_max = ret_deltaNorm_max()
     
  end do all_trajs
  
  !print statistics
  call print_histo()
  write(100,*) '* INPUT******************************************************'
  write(100,*) '* ', N_traj
  write(100,*) '* ', state_init
  write(100,*) '* ', sngl(r_init)
  write(100,*) '* ', sngl(p_init)
  write(100,*) '* ', sngl(mass)
  write(100,*) '* ', sngl(dE), sngl(angle_dp)
  write(100,*) '* ', sngl(upper_limitZ), sngl(lower_limitZ)
  write(100,*) '* ', dt
  write(100,*) '* SUMMARY****************************************************'
  write(100,*) '* deltaE_total_max ', deltaE_max
  write(100,*) '* deltaNorm_total_max ', deltaNorm_max
  write(100,*) '* regular trajectories ', regular_traj
  write(100,*) '* final states', end_state(:)
  if (regular_traj.ne.0) then
     write(100,*) '* final percentage', sngl(dble(end_state(:))/dble(regular_traj))
  else
     write(100,*) '* final percentage', sngl(dble(end_state(:)))
  end if
  write(100,*) '* final population', sngl(population(:)/dble(regular_traj))
  write(100,*) '*************************************************************'
  close(100)

end program surfacehopping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!----------read initial conditions
subroutine read_input_parameters(Ntraj,state_init,r_init,p_init,mass,dE,upper_limitZ,lower_limitZ,dt)
  implicit none

  integer :: Ntraj                                      ! number of trajectories per initial conditions
  integer :: state_init                                 ! initial state
  double precision :: r_init, p_init                    ! initial position and initial momentum
  double precision :: mass                              ! mass of particle
  double precision :: dE, angle_dp                      ! energy width of beam and beam divergence (angle!)
  double precision :: upper_limitZ, lower_limitZ        ! boarders at which trajectory is stopped
  double precision :: dt                                ! time step for integration of cc equations and equation of motion

  open(99,file='input.inp',status='old')
  read(99,*) Ntraj
  read(99,*) state_init
  read(99,*) r_init
  read(99,*) p_init
  read(99,*) mass
  read(99,*) dE
  read(99,*) upper_limitZ, lower_limitZ
  read(99,*) dt
  close(99)

  return
end subroutine read_input_parameters

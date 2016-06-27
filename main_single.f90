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


  ! read input parameters
  call read_input_parameters(N_traj, state_init, r_init, p_init, mass, dE, angle_dp, upper_limitZ, lower_limitZ, dt)

  if (N_traj.ne.1) then
     write(*,*) 'The number of trajectories is larger than 1! Change your input'
     stop
  end if

  call read_input

  call set_nran

  write(*,*)
  write(*,*) "Read input data"
  write(*,*)

  open(99,file='trajectory.out',status='replace')
  open(98,file='velocities.out',status='replace')
  open(97,file='jumps.out',status='replace')
  open(96,file='time_energies.out',status='replace')
  open(95,file='total_energy.out',status='replace')
  open(94,file='a_mat.out',status='replace')
  open(93,file='g_mat.out',status='replace')
  open(92,file='V_mat.out',status='replace')
  open(91,file='b_mat.out',status='replace')
  open(90,file='P_mat.out',status='replace')

  do trajectory = 1,N_traj


     call set_initial_cond(state_init, r_init, p_init, mass, dE, dt)
     call printtraj
     call printstat
     
     do
        call timestep
           call printtraj
           call printstat
           call print_Pmat
!           call printstuff
!        end if

        !projectile is reflected from the surface
        if (ret_actual_posZ().gt.upper_limitZ) exit
        !projectile enters the solid
        if (ret_actual_posZ().le.lower_limitZ) exit

     end do
     
     call print_traj_info()

!     call print_actual_state
!     call print_norm
    
     close(99)
     close(98)
     close(97)
     close(96)
     close(95)
     close(94)
     close(93)
     close(92)
     close(91)
     close(90)
  end do
     
   end program surfacehopping
   


!----------read initial conditions
subroutine read_input_parameters(Ntraj,state_init,r_init,p_init,mass,dE,angle_dp,upper_limitZ,lower_limitZ,dt)
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

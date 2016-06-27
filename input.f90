module input
implicit none
save

public :: N_states, read_input, total_energy, adiabatic_coupling

private
integer, parameter :: N_input=212 ! number of data points for total energies and the nonadiabatic coupling vector
integer, parameter :: N_states=2  ! number of states included in calculation
integer, parameter :: N_coupl=1   ! number of all possible combination of states that can couple

double precision :: R_grid(N_input) ! R grid for total energies and nonadiabatic coupling vector
double precision :: tot_energy(N_input,N_states) ! total energy curves for all included states
double precision :: nad_coupling(N_input,N_coupl) ! nonadiabatic coupling for all possible combinations of states
double precision :: tot_eNEW(N_input,N_states),nad_cNEW(N_input,N_coupl) ! help arrays for preprocessing the input data with spline

contains

!---------- reads input files for total energies and the nonadiabatic coupling vector
  subroutine read_input
    implicit none

    double precision :: dydR1,dydR2
    integer :: i,j

    !read in total energie curves and nonadiabatic coupling
    open(99,file='Vadiabatic.inp',status='old')
!    open(98,file='Coupladiabatic.inp',status='old')

    do i=1,N_input
       read(99,*) R_grid(i),tot_energy(i,:)
!       read(98,*) R_grid(i),nad_coupling(i,:)
    end do

    close(99)
!    close(98)

    !preprocess the total energy curves
    do i=1,N_states
       dydR1=(tot_energy(2,i)-tot_energy(1,i))/(R_grid(2)-R_grid(1))
       dydR2=(tot_energy(N_input,i)-tot_energy(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
       call spline(R_grid,tot_energy(:,i),N_input,dydR1,dydR2,tot_eNEW(:,i))
    end do
!for the coupling a fit formular is used (see function adiabatic_coupling)

!    !preprocess the nonadiabatic coupling
!    do i=1,N_coupl
!       dydR1=(nad_coupling(2,i)-nad_coupling(1,i))/(R_grid(2)-R_grid(1))
!       dydR2=(nad_coupling(N_input,i)-nad_coupling(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
!       call spline(R_grid,nad_coupling(:,i),N_input,dydR1,dydR2,nad_cNew(:,i))
!    end do

  end subroutine read_input


!---------- interpolates the total energy curve for the state "state"
  function total_energy(r,state)
    implicit none
    
    integer :: state
    double precision :: total_energy,help,r

    if (r.lt.R_grid(1)) then
       print*, "WARNING: r is lower than smallest point in R_grid. r=",r
       stop
    end if

    if (r.gt.R_grid(N_input)) then
       total_energy=tot_energy(N_input,state)
       return
    end if

    call splint(R_grid,tot_energy(:,state),tot_eNEW(:,state),N_input,r,help)
    total_energy=help

  end function total_energy


!---------- interpolates the nonadiabatic coupling for the state "state"
  function adiabatic_coupling(r,state)
    implicit none

    double precision :: adiabatic_coupling,r
    !fit parameters of lorenzian  
    double precision :: x0,g,c,pi
    integer :: state
    
    x0 = 6.17124d0
    g = 0.198213d0
    c = 1.51499d0
    pi = acos(-1.d0)
    
    adiabatic_coupling = c/(pi*g*(1+((r-x0)/g)**2))

!original function
!    integer :: state
!    double precision :: adiabatic_coupling,help,r
!
!    if (r.lt.R_grid(1)) then
!       print*, "WARNING: r is lower than smallest point in R_grid"
!       stop
!    end if
!
!    if (r.gt.R_grid(N_input)) then
!       adiabatic_coupling=nad_coupling(N_input,state)
!       return
!    end if
!
!    call splint(R_grid,nad_coupling(:,state),nad_cNEW(:,state),N_input,r,help)
!    adiabatic_coupling=help

    

  end function adiabatic_coupling


!---------- create test output for interpolation of the total energies and the nonadiabatic coupling terms
!  subroutine create_testout
!    implicit none
!
!    integer :: j
!    double precision :: xindex,yindex
!
!    open(77,file='testout_energy.out',status='replace')
!    open(88,file='testout_nadcoupling.out',status='replace')
!
!125 format (5(F15.10,X))
!
!    do xindex=0.6d0,10.d0,0.1d0
!       
!       write(77,*) xindex,(total_energy(xindex,j),j=1,N_states)
!       write(88,125) xindex,(adiabatic_coupling(xindex,j),j=1,N_coupl)
!
!    end do
!
!    close(77)
!    close(88)
!
!  end subroutine create_testout


!---------- preprocesses the input data for interpolation with subroutine splint
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=50000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      end SUBROUTINE spline


!---------- interpolation routine
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) then
        write(*,*) 'bad xa input in splint, program will stop'
        stop
     end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end SUBROUTINE splint

end module input

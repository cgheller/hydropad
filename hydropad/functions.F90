#include "hydrompi.def"

SUBROUTINE timedet
USE dimension
USE scalar
USE matrix
USE mpi_inc
IMPLICIT NONE
real*8 rdtc

at=1.0
ath=1.0
dat=0.0
told = t
dtold = dt
!
! Courant condition
!
write(*,*)"velmax = ", velmax
rdtc=velmax/(cour*dx)
if(velmax .EQ. 0.0)then
  dt=0.1
else
  dt = 1.0/rdtc
endif
t = t+dt
dth=0.5*dt
thalf=told+dth
rdtath=dt/ath


END SUBROUTINE timedet


SUBROUTINE timedet_old
!
USE dimension
USE scalar
USE matrix
USE mpi_inc
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
REAL*8, PARAMETER :: dt0=0.000000001
REAL*8 :: rdt,rdtc,dta
REAL*8 :: rdtrho
REAL*8 :: x,xxx
	REAL*8 :: vpart_max_pe,vpart_max
	REAL*8 :: derivs
!
! calculate quantities at t^n
!
told=t
dtold=dt
if(nstep.gt.1)then
   rdt=5.0/(6.0*dt)
else
   rdt=0.0
endif
!
#ifdef COSMO
if(omega_m.eq.1.0) then
  at=t**(2.0/3.0)
  dat=(2.0/3.0)*t**(-1.0/3.0)
  dta=((1.0+ca)*at)**(3.0/2.0)-t
else
  dat=datnew
  at=atnew
  dta=ca*at/dat
end if
!
! find maximum particles velocity
!
vpart_max_pe=maxval(v1**2+v2**2+v3**2)
vpart_max_pe=sqrt(vpart_max_pe)
#ifdef USEMPI
        CALL MPI_Allreduce(vpart_max_pe,vpart_max,&
                          1,MPI_DOUBLE_PRECISION,&
                          MPI_MAX,MPI_COMM_WORLD,ierr)
#else
vpart_max=vpart_max_pe
#endif
#else
dta=1e10
vpart_max=0.0
at=1.0
#endif
!
! calculate new timestep
!
velmax=max(velmax,vpart_max)
rat=1.0/at
rdtc=velmax/(cour*dx*at)
dt=1.0/rdtc
!CLA	dt=1.0/max(rdt,rdtc)
dt=min(dt,dta)
if(dt .ge. 10.0)dt=0.1
!CLA reintroduce the following I/O
#ifdef CLA
if(nstep.gt.1)then
  if(dt.lt.dt0)then
     CALL outdumpbm
     CALL outdumpdm
!
#ifdef USEMPI
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
    if(mype.eq.0)then
      open(unit=10,file='terror.dat')
      write(10,*)'dt is less than...',dt0
      write(10,*)'Simulation is dumped'
      write(10,*)'dt = ',dt
      write(10,*)'dt_old = ',dtold
      close(10)
      stop
    endif
  endif
endif
#endif
!!!	dt=min(dt,dt0)
if(dt.eq.dt0)write(*,*)'dt = dt0'
t=t+dt
write(*,*)"---------> ", told, t, dt
if(t.gt.tstop)then
   t=tstop
   dt=t-told
endif
dth=0.5*dt
thalf=told+dth
!
! calculate quantities at t^{n+1} and t^{n+1/2}
!
#ifdef COSMO
if(omega_m.eq.1.0) then
  atnew=t**(2.0/3.0)
  ath=thalf**(2.0/3.0)
  dath=(2.0/3.0)*thalf**(-1.0/3.0)
  datnew=(2.0/3.0)*t**(-1.0/3.0)
else
  if(mype.eq.0)then
    CALL rk4(at,dat,1,told,dt,atnew)
    CALL rk4(at,dat,1,told,dth,ath)
  endif
!
#ifdef USEMPI
CALL MPI_Bcast(atnew,1,MPI_DOUBLE_PRECISION,0,&
               MPI_COMM_WORLD,ierr)
CALL MPI_Bcast(ath,1,MPI_DOUBLE_PRECISION,0,&
               MPI_COMM_WORLD,ierr)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
  datnew=derivs(t,atnew)
  dath=derivs(thalf,ath)
!
endif
! 
rdtath=dt/ath
!
! WARNING: omega is related only to the matter component (Omega_lambda 
! is constant !!!)
!
omega=real(1.0/(1.0+atnew*(1.0-omega_m)/omega_m))
hubble=at/dat
redshiftold=1.0/at-1.0
redshift=1.0/atnew-1.0
#endif
!
if(mype.eq.0)then
  write(1150,1101)t,atnew,redshift,omega,datnew,t**0.666667,&
                  0.6666667*t**(-0.3333333)
1101  format(7(1x,e11.5))
endif
!
END SUBROUTINE timedet_old
!
!********************************************************************
!********************************************************************
!
SUBROUTINE start
!
USE dimension
USE scalar
IMPLICIT NONE
INTEGER :: i,j,k
!
! local variables
!
real*8 dk
!
dk=twopi/xmax
do k=0,n1/2-1
   ak(k)=float(k)
enddo
do k=n1/2,n11
   ak(k)=float(k-n1)
enddo
do k=0,n11
   ak(k)=ak(k)*dk
   akq(k)=ak(k)*ak(k)
enddo

END SUBROUTINE start
!
!********************************************************************
!********************************************************************
!
REAL*8 FUNCTION derivs(timeaux,ataux)
!
USE dimension
USE scalar
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8, INTENT(IN) :: timeaux,ataux
!
derivs=t0h0*sqrt(1+omega_m*(1.0/ataux-1.0)+&
                 omega_lambda*(ataux**2-1.0))

END FUNCTION derivs

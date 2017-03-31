SUBROUTINE initialize_parameters
!
USE dimension
USE matrix
USE vector
USE scalar
USE io_mod
#ifdef TRACERS
USE tracers_mod
#endif
!
USE mpi_inc
!
IMPLICIT NONE
!
! define local variables
!
INTEGER :: i,j,k
REAL*8, EXTERNAL :: s11abf
REAL*8 :: derivs
!
amass=1.0
dx=1.0/real(max(ngridx, ngridy, ngridz))
rdx=1.0/dx
!
! calculate fundamental density parameters:
!
#ifndef COSMO
omega_dm = 0.0
omega_lambda = 0.0
omega_bm = 1.0
boxMpc_over_h = 1.0
hnow = 1.0
#endif
omega_m=omega_dm+omega_bm
omega_tot=omega_dm+omega_bm+omega_lambda
omega_a=omega_m-0.3*omega_tot+0.3
omega0=omega_tot
if(omega_dm.ne.0.0)then
   dmfact=1.0/omega_dm
else
   dmfact=1.0
endif
rdmfact=1.0/dmfact
amass=rdmfact/float(npart)
bmfact=1.0/omega_bm
rbmfact=1.0/bmfact
box=boxMpc_over_h/hnow
dxmpc=box/float(ngridx)
if(omega_m.eq.1.0)then
  t0h0=0.66666666667
else
  xsinh=dble(sqrt(abs(1.0-omega_a)/omega_a))
  if(mype.eq.0)t0h0=0.666667*s11abf(xsinh,ifail)/sqrt(abs(1.0-omega_a))
endif
#ifdef USEMPI
CALL MPI_Bcast(t0h0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

!
! calculate derived physics quantities 
!
r=boltz/mh
rgamma=1.0/gamma
gamma1=(gamma+1.0)/(gamma-1.0)
rgamma1=1.0/(gamma-1)
m=(gamma-1.0)/(2.0*gamma)
rm=1.0/m
gf=(gamma+1.0)/(2.0*gamma)
!
! set cosmology related physical quantitites
!
! a=1 box resolution in cm 
!
xnow=3.086e24*box/float(ngridx)
!
! hubble constant in 1/sec
!
h_sec=hnow/3.0856e17
!
rhonow=3.0*h_sec**2/(25.13274*6.6726e-8)
   tnow=t0h0/h_sec
#ifdef COSMO
if(omega_m.eq.1.0)then
   tstop=(tstop_aux+1.0)**(-1.5)
!
else
   if(tstop.ne.0.0)then
      write(*,*)'final redshift not equal to 0'
      write(*,*)'is not supported for this model'
      stop
   else
      tstop=1.0
   endif
!
endif
#endif
!
! gravitational constant in code units
!
gconst=4.0*pi*(6.6726e-8)*rhonow*tnow*tnow
!
! constant to get physical pressure from code pressure
!
pfact=r*(tnow/xnow)**2*tempnow
!
! constant to get physical temperature (K) from code temperature
!
tfact=(xnow/tnow)**2/r
!
! constant to get physical velocity (cm/sec) from code velocity
!
vnor=1e-5*xnow/tnow
hfact=hfrac**2*rhonow/mh**2
hfact=(hfact*tnow*tnow)*tnow
hfact=(hfact/xnow)/xnow
cmfact=(tnow*tnow/mh)*tnow*hfrac
cmfact=(cmfact/xnow)/xnow
!
! speed of light in code units
!
light=(tnow/xnow)*2.9979e10
!
! useful constant
!
tsux=(tnow/xnow)**2
!
! if not startdump set initial time variables and normalize velocities
!
!
! initialize time parameters
!
#ifdef COSMO
!
if(startdump.ne.1)then
! approximate initial time by a flat model
!
          tin=at**1.5
!
          told=tin
          dt=0.010*tin
          dth=0.5*dt
          dtold=dt
          t=tin+dt
endif
        thalf=told+dt/2.0
!
if(omega_m.eq.1.0) then
           at=told**(2.0/3.0)
           atnew=t**(2.0/3.0)
           dat=(2.0/3.0)*told**(-1.0/3.0)
           ath=(told+dt/2.)**(2.0/3.0)
           dath=(2.0/3.0)*thalf**(-1.0/3.0)
           datnew=(2.0/3.0)*t**(-1.0/3.0)
else
           dat=derivs(told,at)
           if(mype.eq.0)then
             CALL rk4(at,dat,1,told,dth,ath)
             CALL rk4(at,dat,1,told,dt,atnew)
           endif
!
#ifdef USEMPI
           CALL MPI_Bcast(ath,1,MPI_DOUBLE_PRECISION,0,&
                          MPI_COMM_WORLD,ierr)
           CALL MPI_Bcast(atnew,1,MPI_DOUBLE_PRECISION,0,&
                          MPI_COMM_WORLD,ierr)
#endif
           dath=derivs(thalf,ath)
           datnew=derivs(t,atnew)
!
endif
!
rat=1.0/at
dtat=dt*rat
rdtath=dt/ath
omega=real(1.0/(1.0+atnew*(1.0-omega_m)/omega_m))
hubble=at/dat
redshiftold=1.0/at-1.0
redshift=1.0/atnew-1.0
!
#else
!
tin=0.0
told=tin
dt=dtinit
dth=0.5*dt
dtold=dt
t=tin+dt
thalf=told+dt/2.0
at=1.0
atnew=at
dat=0.0
dtat=dt
ath=1.0
dath=0.0
datnew=0.0
rat=1.0/at
rdtath=dt/ath
omega=1.0
hubble=1.0
redshiftold=0.0
redshift=0.0
thalf=told+dth
!
#endif

END SUBROUTINE initialize_parameters

#include "hydrompi.def"
!
SUBROUTINE initialize(startdump)
!
	USE dimension
	USE matrix
	USE vector
	USE scalar
	USE mpi_inc
	USE io_mod
	USE tracers_mod
!
	IMPLICIT NONE
!
! local variables
!
	INTEGER, INTENT(IN) :: startdump
	INTEGER :: n_tot
	INTEGER :: iunit,iunit1
	INTEGER :: iloop,nloop
	REAL*8 :: dtaux,vfact
	REAL*8 :: dum1(1),dum2(1),dum3(1)
	REAL*8, EXTERNAL :: s11abf
	REAL*8 :: derivs

!
	dum1=0.0
	dum2=0.0
	dum3=0.0
	n_tot=1
	initboxyn=0
	initcheckyn=0
	initproject=0
!
	if(omega_m.eq.1.0)then
	  t0h0=0.66666666667
	else
          xsinh=dble(sqrt(abs(1.0-omega_a)/omega_a))
          if(mype.eq.0)then
            t0h0=0.666667*s11abf(xsinh,ifail)/sqrt(abs(1.0-omega_a))
          endif
#ifdef USEMPI
          CALL MPI_Bcast(t0h0,1,MPI_DOUBLE_PRECISION,0,&
                         MPI_COMM_WORLD,ierr)
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	endif
!
! initialize system if starts from dumped initial conditions
!
	if(startdump.eq.1)then
	   CALL indumpbm
	   dth=0.5*dt
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	   if(npes.gt.1)then
	     CALL initgold
!
#ifdef USEMPI
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	   endif
!
	   CALL indumpdm
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
!
! initialize system if start from initial conditions
!
	else
!
! calculate initial velocities and postions of particles
!
	   noutput=0
	   CALL initialdm
!
#ifdef USEMPI
           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	endif
!
! calculate normalization factors >
! a=1 box size in cm 
!
	xnow=3.086e24*box	
!
! hubble constant in 1/sec
!
	h_sec=hnow/3.0856e17
!
        rhonow=3.0*h_sec**2/(25.13274*6.6726e-8)
!
        if(omega_m.eq.1.0)then
             tnow=t0h0/h_sec
             tstop=(tstop_aux+1.0)**(-1.5)
!
        else
	     tnow=t0h0/h_sec
!
	     if(tstop.ne.0.0)then
	       write(*,*)'final redshift not equal to 0'
	       write(*,*)'is not supported for this model'
	       stop
	     else
	       tstop=1.0
	     endif
!
	endif
!
! gravitational constant in code units
!
	   gconst=4.0*pi*(6.6726e-8)*rhonow*tnow*tnow 
!
! size in cm of the cell side
!
	   xnow=xnow/float(ngrid)
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
	if(startdump.ne.1)then
!
! initialize time parameters
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
	   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	   dath=derivs(thalf,ath)
	   datnew=derivs(t,atnew)
!
        end if
!
	rat=1.0/at
	rdtath=dt/ath
!
! calculate initial dark matter density field
!
	if(startdump.ne.1)then
!
! normalize velocities from km/sec to code units
!
           vfact=rat*1e5*tnow/xnow
    
           write(*,*)"TNOW/XNOW = ",tnow/xnow
           v1=vfact*v1
           v2=vfact*v2
           v3=vfact*v3
!
           CALL dynamics
!
#ifdef USEMPI
           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! calculate initial baryonic matter density and pressure fields
!
	   CALL initialbm
!
#ifdef USEMPI
           CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	endif
!
! 
        omega=real(1.0/(1.0+atnew*(1.0-omega_m)/omega_m))
	hubble=at/dat
        redshiftold=1.0/at-1.0
        redshift=1.0/atnew-1.0
!
! initialize tracers
!
	checkstarttrack=0
!	if(redshiftold.le.zstarttrack)CALL tracers_init(startdump)
!
END SUBROUTINE initialize

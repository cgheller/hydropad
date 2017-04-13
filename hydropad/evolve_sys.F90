SUBROUTINE evolve_sys
!
! evolve the system from t^n to t^n+1
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
USE nbody_mod
#ifdef PROFILE
USE profile
#endif
#ifdef TRACERS
USE tracers_mod
#endif
!
! local variables
!
REAL(KIND=8) :: tm1,tm2,tm3,pmax1,pmax2,pmax3
INTEGER :: i,j,k
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: statcomm
INTEGER, ALLOCATABLE, DIMENSION(:) :: request
INTEGER :: nmoved,nreq
!
nreq = 4*npes
ALLOCATE(request(nreq))
ALLOCATE(statcomm(MPI_STATUS_SIZE,nreq))
!
nrot=mod(nstep,6)+1
if(nstep == 1)nrot=1
!
! calculate the DM density field, the gravitational potential and forces
! at the beginning of the step - move dark matter particles
!
#ifdef GRAVITY
!
! in dynamics right now ther is the calculation of density, of gravitational
! field and of the nbody. Furthermore there is all the communication for the particles. 
! Must be split!!!
!
#endif
!
! exchange ghost regions
!
!CALL exchange_mesh
!
! call hydro solver
!
#ifdef HYDRO
#ifndef STENCIL
CALL savetnp1
#endif
CALL ppm_lxlylz
#endif
!
! save t^n values of hydro variables
!
!CALL savetn
!
! Nbody
!
#ifdef NBODY
#ifdef USEMPI
CALL exchange_parts(npartpe,ppos,pvel,request,nreq)
nmoved=nrecv
#endif
CALL kick(npartpe,ppos,pvel,nx,ny,nz,gforce)
#ifdef USEMPI
CALL MPI_Waitall(nreq,request,statcomm,ierr)
CALL kick(nmoved,xrecv,vrecv,nx,ny,nz,gforce)
#endif
CALL drift(npartpe,ppos,pvel)
#ifdef USEMPI
CALL drift(nmoved,xrecv,vrecv)
#endif
! here we have to calculate the density and the gravitational field
CALL kick(npartpe,ppos,pvel,nx,ny,nz,gforce)
#ifdef USEMPI
CALL kick(nmoved,xrecv,vrecv,nx,ny,nz,gforce)
CALL exchangeback_parts(npartpe,ppos,pvel)
CALL dealloc_particles
#endif
#endif
!
! expansion step
!
#ifdef COSMO
CALL expand
#endif
!
! calculate average temperature and max pressure
!
!!!!!!!!CALL calculate_thermo(tm1,pmax1)
!
!
! cooling
!
#ifdef COOLING
if(ncool.eq.1)then
  CALL tablecool
  CALL calculate_thermo(tm2,pmax2)
endif
#endif
!
! external supernovae heating
!
#ifdef SUPERNOVAE
if(nuv.eq.1)CALL heating
CALL calculate_thermo(tm3,pmax3)
#endif
#ifdef TRACERS
if(redshiftold.le.zstarttrack)then
  if(checkstarttrack.eq.0)CALL tracers_init(startdump)
  CALL tracers
  if(mod(nstep,troutstep).eq.0)CALL outtracers(0)
endif
#endif
!
! calculate maximum velocity
!
CALL speed
!
! calculate the new timestep
!
#ifdef COSMO
CALL timedet_cosmo
#else
CALL timedet
#endif
!
! print mean quantities
!
!!!!!!CALL means(tm1,tm2,tm3,pmax1,pmax2,pmax3)
!
DEALLOCATE(request)
DEALLOCATE(statcomm)
!
END SUBROUTINE evolve_sys

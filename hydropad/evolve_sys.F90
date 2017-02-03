#include "hydrompi.def"
!
SUBROUTINE evolve_sys(startdump)
!
! evolve the system from t^n to t^n+1
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
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
REAL(KIND=8) ::  ti,tf
INTEGER :: startdump,dump
INTEGER :: i,j,k
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
call dynamics
#endif
!
#ifdef COSMO
!
! CLA: THIS HAS TO BE MOVED IN init_sys 
!
!!if(nstep.eq.1)then
!!   call initvel
!!   nrot=1
!!endif
#endif
!
! exchange ghost regions
!
!call exchange_mesh
!
! search for shocked regions
! CLA: MUST BE REWRITTEN
!
call shsearch
!
! call hydro solver
!
call savetnp1
call ppm_lxlylz
!
! save t^n values of hydro variables
!
!call savetn
!
! expansion step
!
#ifdef COSMO
call expand
#endif
!
! calculate average temperature and max pressure
!
call calculate_thermo(tm1,pmax1)
!
nrot=mod(nstep,6)+1
!
! cooling
!
#ifdef COOLING
if(ncool.eq.1)then
  call tablecool
  call calculate_thermo(tm2,pmax2)
endif
#endif
!
! external supernovae heating
!
#ifdef SUPERNOVAE
if(nuv.eq.1)call heating
call calculate_thermo(tm3,pmax3)
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
call speed
!
! calculate the new timestep
!
call timedet
!
! print mean quantities
!
call means(tm1,tm2,tm3,pmax1,pmax2,pmax3)
!
END SUBROUTINE evolve_sys

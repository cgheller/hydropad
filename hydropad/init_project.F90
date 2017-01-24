#include "hydrompi.def"
!
SUBROUTINE init_project
!
        USE dimension
        USE scalar
        USE mpi_inc
	USE xproj
	USE io_mod

	IMPLICIT NONE
	INTEGER :: nchan1,ntemp1
	INTEGER :: flag_telescope
!
	initproject=1
	open(unit=12,file='project.in')
	read(12,*)nprj
	read(12,*)zstart

	READ(12,*) powmin
  	READ(12,*) powmax
  	READ(12,*) nchan1
  	READ(12,*) de
  	READ(12,*) Tmin
  	READ(12,*) dTT
  	READ(12,*) ntemp1
  	READ(12,*) zmet

	read(12,'(a)')filepow
	read(12,*)flag_telescope

	close(12)
!
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	Tmax=10**(log10(Tmin)+dTT*(ntemp-1.))
!
	call loadpow
!
	call get_limits(flag_telescope,ibin1,ibin2)
!	
END SUBROUTINE init_project

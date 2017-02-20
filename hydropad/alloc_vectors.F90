SUBROUTINE alloc_vectors(dgrid)
!
USE dimension
USE vector
USE mpi_inc
!
IMPLICIT NONE
INTEGER :: error
INTEGER :: dgrid

ALLOCATE (cho(dgrid),STAT=error)
ALLOCATE (nes(dgrid),STAT=error)
ALLOCATE (pres(dgrid),STAT=error)
ALLOCATE (rho(dgrid),STAT=error)
ALLOCATE (vx(dgrid),STAT=error)
ALLOCATE (vy(dgrid),STAT=error)
ALLOCATE (vz(dgrid),STAT=error)
ALLOCATE (ghalf(dgrid),STAT=error)
ALLOCATE (g(dgrid),STAT=error)
ALLOCATE (c(dgrid),STAT=error)
ALLOCATE (eint(dgrid),STAT=error)
ALLOCATE (etot(dgrid),STAT=error)
ALLOCATE (rhoold1d(dgrid),STAT=error)
ALLOCATE (vxold1d(dgrid),STAT=error)
ALLOCATE (phi1(dgrid),STAT=error)
ALLOCATE (phi0(dgrid),STAT=error)
ALLOCATE (flt(dgrid),STAT=error)
ALLOCATE (nu(dgrid),STAT=error)
ALLOCATE (sk2(dgrid),STAT=error)
!
!$acc kernels
cho=0.0
nes=0.0
pres=0.0
rho=0.0
vx=0.0
vy=0.0
vz=0.0
ghalf=0.0
g=0.0
c=0.0
eint=0.0
etot=0.0
rhoold1d=0.0
vxold1d=0.0
phi1=0.0
phi0=0.0
flt=0.0
nu=0.0
sk2=0.0
!$acc end kernels
!
END SUBROUTINE alloc_vectors

MODULE vector
!
!USE dimension
IMPLICIT NONE
SAVE
!
! define basic one dimensional vectors
!
REAL*8, DIMENSION(:), ALLOCATABLE :: cho
REAL*8, DIMENSION(:), ALLOCATABLE :: nes
!
REAL*8, DIMENSION(:), ALLOCATABLE :: pres,rho,vx,vy,ghalf
REAL*8, DIMENSION(:), ALLOCATABLE :: vz,g,c,eint,etot
REAL*8, DIMENSION(:), ALLOCATABLE :: rhoold1d,vxold1d,phi1,phi0
REAL*8, DIMENSION(:), ALLOCATABLE :: flt,nu,sk2
!
END MODULE vector

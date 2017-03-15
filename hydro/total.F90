REAL*8 FUNCTION total(pres,rho,vx,vy,vz)
!
! calculate internal energy density
!
!USE dimension
!USE vector
USE scalar
!$acc routine seq
!
! local variables
!
IMPLICIT NONE
REAL(KIND=8) :: pres,rho,vx,vy,vz
!
total=pres*rgamma1+0.5*rho*(vx*vx+vy*vy+vz*vz)
!
END FUNCTION total

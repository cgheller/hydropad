SUBROUTINE flatten5(abar,flatyn,ncell)

IMPLICIT NONE
!$acc routine seq
REAL(KIND=8) :: abar(ncell)
INTEGER :: flatyn
REAL(KIND=8) :: gradleft,gradright,gradmax,mingradflat
INTEGER :: i,ncell,icenter
!
flatyn = 0
!
gradmax = 0.0
do i = 1,4
  icenter = i+1
  gradleft   = abs((abar(icenter)-abar(icenter-1))/min(abar(icenter),abar(icenter-1)))
  gradmax = max(gradmax,gradleft)
enddo
!
if(gradmax >= mingradflat)flatyn = 1

END SUBROUTINE flatten5

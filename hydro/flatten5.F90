SUBROUTINE flatten5(abar,flatyn,ncell)
!$acc routine seq

USE scalar
IMPLICIT NONE
REAL(KIND=8) :: abar(ncell)
INTEGER :: flatyn
REAL(KIND=8) :: gradleft,gradright,gradmax
INTEGER :: i,ncell,icenter
!!!!$acc update device (mingradflat)
!!!!$acc data present (abar)
!
gradmax = 0.0
i=1
!do i = 1,4
  icenter = i+1
  if(abar(icenter)*abar(icenter-1) .NE. 0.0)&
         gradleft = abs((abar(icenter)-abar(icenter-1))/min(abar(icenter),abar(icenter-1)))
  gradmax = max(gradmax,gradleft)
!enddo
i=2
  icenter = i+1
  if(abar(icenter)*abar(icenter-1) .NE. 0.0)&
         gradleft = abs((abar(icenter)-abar(icenter-1))/min(abar(icenter),abar(icenter-1)))
  gradmax = max(gradmax,gradleft)
i=3
  icenter = i+1
  if(abar(icenter)*abar(icenter-1) .NE. 0.0)&
         gradleft = abs((abar(icenter)-abar(icenter-1))/min(abar(icenter),abar(icenter-1)))
  gradmax = max(gradmax,gradleft)
i=4
  icenter = i+1
  if(abar(icenter)*abar(icenter-1) .NE. 0.0)&
         gradleft = abs((abar(icenter)-abar(icenter-1))/min(abar(icenter),abar(icenter-1)))
  gradmax = max(gradmax,gradleft)
!
if(gradmax >= mingradflat)flatyn = 1
!!!!!$acc end data

END SUBROUTINE flatten5

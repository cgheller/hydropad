SUBROUTINE means(tm1,tm2,tm3,pmax1,pmax2,pmax3)
!
! store of some mean quantities
!
USE dimension
USE vector
USE scalar
!
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
REAL*8, INTENT(OUT) :: tm1,tm2,tm3,pmax1,pmax2,pmax3
REAL*8 :: sqsigmatot,sqsigma,sqsigma1,dddmax,dddmax1
!
common/out/sqsigma,sqsigma1,dddmax,dddmax1
!
if(mype.eq.0)then
open(900,file='av_density.dat',position='append')
write(900,10)redshiftold,sqsigma,sqsigma1,dddmax,dddmax1
close(900)
open(1000,file='av_temp.dat',position='append')
write(1000,11)redshiftold,tm1,tm2,tm3,pmax1,pmax2,pmax3
close(1000)
endif
10  format(5(1x,e11.5))
11  format(7(1x,e11.5))
!
END SUBROUTINE means

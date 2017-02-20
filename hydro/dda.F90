SUBROUTINE dda(a,eta,ajld,ajrd,dmalf,dmarg)
!
USE dimension
USE scalar
!$acc routine seq
!
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8 :: a(5)
REAL*8 :: eta,ajld,ajrd,dmalf,dmarg,deltam2
REAL*8 :: deltap2,agrad,mmmax,mmmin
!
! linear correction
!
mmmax=0.0
mmmin=1.0
deltam2=0.16666667*(a(3)-2.0*a(2)+a(1))/(dx*dx)
deltap2=0.16666667*(a(5)-2.0*a(4)+a(3))/(dx*dx)
agrad=abs(a(4)-a(2))-0.010*min(abs(a(4)),abs(a(2)))
if(-deltam2*deltap2.gt.0.0.and.agrad.gt.0.0)then
    eta=-(deltap2-deltam2)*dx*dx/(a(4)-a(2))
    eta=max(mmmax,min(20.0*(eta-0.050),mmmin))
else
    eta=0.0
endif
!
! left linear interpolation 
!
ajld=a(2)+0.50*dmalf
!
! right linear interpolation
!
ajrd=a(4)-0.50*dmarg
!
END SUBROUTINE dda

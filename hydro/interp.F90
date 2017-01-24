SUBROUTINE interp(a,ajl,ajr,deltaa1,a6,flat,ndda,addp,addrho)
!
! this routine evaluates the interpolation function of the hydro
! quantities
!
USE scalar
!
IMPLICIT NONE
INTEGER ndda
INTEGER :: i,j,k
REAL*8 :: a(5)
REAL*8 :: ajl(2),ajr(2),deltaa1(2),a6(2)
REAL*8 :: flat,ajld,ajrd,eta,addrho,addp,ajside
REAL*8 :: dmaj,dmajp1,dmajp2
REAL*8 vi,vip1,atl,atr,acl,acr 
!
!
! compute left and right limits of the hydro functions for the
! (i+1)-zone
!
dmaj=ajside(a(1),a(2),a(3))
dmajp1=ajside(a(2),a(3),a(4))
dmajp2=ajside(a(3),a(4),a(5))
ajl(2)=a(2)+0.50*(a(3)-a(2))+(dmaj-dmajp1)*0.16666667
ajr(2)=a(3)+0.50*(a(4)-a(3))+(dmajp1-dmajp2)*0.16666667
if(ndda.eq.1.and.addrho.ge.addp)then
   call dda(a,eta,ajld,ajrd,dmaj,dmajp2)
!          eta=0.0
   ajl(2)=(1.0-eta)*ajl(2)+eta*ajld
   ajr(2)=(1.0-eta)*ajr(2)+eta*ajrd
endif
deltaa1(2)=ajr(2)-ajl(2)
a6(2)=6*(a(3)-0.50*(ajl(2)+ajr(2)))
if((ajr(2)-a(3))*(a(3)-ajl(2)).le.0)then
       ajl(2)=a(3)
       ajr(2)=a(3)
else if(deltaa1(2)*a6(2).gt.deltaa1(2)**2)then
       ajl(2)=3.0*a(3)-2.0*ajr(2)
else if(deltaa1(2)*a6(2).lt.(-(deltaa1(2)**2)))then
       ajr(2)=3.0*a(3)-2.0*ajl(2)
endif
!flatttt	flat=0.0
ajl(2)=flat*a(3)+(1.0-flat)*ajl(2)
ajr(2)=flat*a(3)+(1.0-flat)*ajr(2)
deltaa1(2)=ajr(2)-ajl(2)
a6(2)=6*(a(3)-0.50*(ajl(2)+ajr(2)))
!
END SUBROUTINE interp

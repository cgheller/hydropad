SUBROUTINE correction(addp,addrho,a,ajl,ajr,deltaa1,a6,ccc,vi,vip1,atl,atr,acl,acr,flat,ndda)
!
! this routine evaluates the corrected values for the hydro 
! variables that are used to compute the initial left and
! right states of the fluid for the Riemann problem.
!
USE dimension
USE scalar
!
IMPLICIT NONE
INTEGER :: ndda
INTEGER :: i,j,k
REAL*8  :: a(5)
REAL*8  :: acl(3),acr(3)
REAL*8  :: ajl(2),ajr(2),deltaa1(2),a6(2),ccc(2)
REAL*8  :: atl,atr,fl,fr,y,vi,vip1,cv
REAL*8  :: flat,ajld,ajrd,eta,addrho,addp,ajside
REAL*8  :: dmaj,dmajp1,dmajp2,mmmax

!
! compute left and right limits of the hydro functions for the
! (i+1)-zone
!
mmmax=0.0
dmaj=ajside(a(1),a(2),a(3))
dmajp1=ajside(a(2),a(3),a(4))
dmajp2=ajside(a(3),a(4),a(5))
ajl(2)=a(2)+0.50*(a(3)-a(2))+(dmaj-dmajp1)*0.16666667
ajr(2)=a(3)+0.50*(a(4)-a(3))+(dmajp1-dmajp2)*0.16666667
if(ndda.eq.1.and.addrho.ge.addp)then
   call dda(a,eta,ajld,ajrd,dmaj,dmajp2)
   !eta=0.0
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
ajl(2)=flat*a(3)+(1.0-flat)*ajl(2)
ajr(2)=flat*a(3)+(1.0-flat)*ajr(2)
deltaa1(2)=ajr(2)-ajl(2)
a6(2)=6*(a(3)-0.50*(ajl(2)+ajr(2)))
!
! compute a first approximation of the left and right states
! for the fluid
!
cv=dt*(vi+ccc(1))*rat       
y=max(cv,mmmax)
atl=fl(ajr,deltaa1,a6,y)
cv=-dt*(vip1-ccc(2))*rat 
y=max(cv,mmmax)
atr=fr(ajl,deltaa1,a6,y)
!
! compute the corrections to the previous approximation
! (1=-, 2=0, 3=+)
! left corrections:
!
y=dt*(vi-ccc(1))*rat
acl(1)=fl(ajr,deltaa1,a6,y)
y=dt*vi*rat
acl(2)=fl(ajr,deltaa1,a6,y)
y=dt*(vi+ccc(1))*rat
acl(3)=fl(ajr,deltaa1,a6,y)
!
! right corrections:
!
y=-dt*(vip1-ccc(2))*rat
acr(1)=fr(ajl,deltaa1,a6,y)
y=-dt*vip1*rat
acr(2)=fr(ajl,deltaa1,a6,y)
y=-dt*(vip1+ccc(2))*rat
acr(3)=fr(ajl,deltaa1,a6,y)
!
END SUBROUTINE correction
!
!*************************************************************
!*************************************************************
!
REAL*8 FUNCTION fl(ajr,deltaa1,a6,y)
!
USE dimension
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8 :: ajr(2),deltaa1(2),a6(2)
REAL*8 :: y,x
!
x=y/dx
fl=ajr(1)-0.50*x*(deltaa1(1)-(1.0-2.0*x*0.33333333)*a6(1))
!
END FUNCTION fl
!
!*************************************************************
!*************************************************************
!
REAL*8 FUNCTION fr(ajl,deltaa1,a6,y)
!
USE dimension
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8 :: ajl(2),deltaa1(2),a6(2)
REAL*8 :: y,x
!
x=y/dx
fr=ajl(2)+0.50*x*(deltaa1(2)+(1.0-2.0*x*0.33333333)*a6(2))
!
END FUNCTION fr
!
!*************************************************************
!*************************************************************
!
REAL*8 FUNCTION beta(ct,ut,uc,pt,pc,gc,sign)
!
USE dimension
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8 :: ct,ut,uc,pt,pc,gc,sign,gvc
!
gvc=(-gc+dat*uc)*rat
beta=-sign*( (ut-uc)+sign*&
     (pt-pc)/ct+sign*dt*(2.0*dat*pc/(ct*at))+&
     dt*gvc)/(2*ct)
!
END FUNCTION beta
!
!*************************************************************
!*************************************************************
!
REAL*8 FUNCTION beta0(ct,rhot,rhoc,pt,pc)
!
USE dimension
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8 :: ct,rhot,rhoc,pt,pc
!
beta0=(pt-pc)/(ct*ct)+1.0/rhot-1.0/rhoc+&
      (2.0*dat*pc)/(at*ct*ct)*dt
!
END FUNCTION beta0
!
!*********************************************************************
!*********************************************************************
!
REAL*8 FUNCTION ajside(am1,a,ap1)
!
IMPLICIT NONE
INTEGER :: i,j,k
REAL*8 :: am1,a,ap1,da
!
if((ap1-a)*(a-am1).gt.0.0)then
  da=0.50*(ap1-am1)
  ajside=min(abs(da),2.0*abs(a-am1),2.0*abs(ap1-a))*da/abs(da)
else
  ajside=0.0
endif
!
END FUNCTION ajside

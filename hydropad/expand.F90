SUBROUTINE expand
!
! compute the contribution of the expansion on the hydro variables
!
USE dimension
USE matrix
USE vector
USE scalar
!
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
INTEGER :: nnncheck
REAL*8  :: afact,afact1,afact2,ppp
!
afact=dt*(dath/ath)
afact1=1.0/(1.0+afact)
afact2=1.0/(2.0+afact)
!
do k=1,nz
 do j=1,ny
  do i=1,nx
    vx3d(i,j,k)=(2.0*vx3d(i,j,k)-afact*vxold(i,j,k))*afact2
    vy3d(i,j,k)=(2.0*vy3d(i,j,k)-afact*vyold(i,j,k))*afact2
    vz3d(i,j,k)=(2.0*vz3d(i,j,k)-afact*vzold(i,j,k))*afact2
    ppp=p3d(i,j,k)
    p3d(i,j,k)=(p3d(i,j,k)-afact*pold3d(i,j,k))*afact1
!
   if(p3d(i,j,k).le.0.0)then
       if(ppp.gt.0.0)then
          p3d(i,j,k)=ppp
       else
          p3d(i,j,k)=pold3d(i,j,k)
       endif
   endif
!
   enddo
  enddo
enddo
!
END SUBROUTINE expand

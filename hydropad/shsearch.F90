#include "hydrompi.def"

SUBROUTINE shsearch
!
! search for strongly compressing regions (possible shocks)
! 
USE dimension
USE matrix
USE scalar
USE mpi_inc
!
! local variables
!
IMPLICIT NONE
!
INTEGER :: indexm,indexp,nslice
INTEGER :: i,j,k
REAL*8  :: pgradx,vgradx,pgrady,vgrady,pgradz,vgradz,vgrad,pgrad
!
pgrad=0.0
!
do k=nbound+1,nz-nbound
  do j=nbound+1,ny-nbound
     do i=nbound+1,nx-nbound
      nes3d(i,j,k)=1.0
!
! compute pressure and velocity gradients in the x-direction:
!
      indexp=i+1
      indexm=i-1
      pgradx=abs(p3d(indexp,j,k)-p3d(indexm,j,k))
      vgradx=vx3d(indexp,j,k)-vx3d(indexm,j,k)
!
! y-direction:
!
      indexp=j+1
      indexm=j-1
      pgrady=abs(p3d(i,indexp,k)-p3d(i,indexm,k))
      vgrady=vy3d(i,indexp,k)-vy3d(i,indexm,k)
!
! z-direction:
!
      indexp=k+1
      indexm=k-1
      pgradz=abs(p3d(i,j,indexp)-p3d(i,j,indexm))
      vgradz=vz3d(i,j,indexp)-vz3d(i,j,indexm)
!
      vgrad=vgradx+vgrady+vgradz
      pgrad=sqrt(pgradx*pgradx+pgrady*pgrady+pgradz*pgradz)/p3d(i,j,k)
      if(vgrad.lt.0.0.and.pgrad.gt.eta1)nes3d(i,j,k)=2.0
    enddo
  enddo
enddo
!CLA
nes3d=1.0
!cho3d=1.0
!
END SUBROUTINE shsearch

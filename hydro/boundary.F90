SUBROUTINE boundary(flat)
!
! compute boundary values of the hydro quantities with periodic
! boundary conditions
!
USE dimension
USE vector
USE scalar
!
! local variables
!
INTEGER :: i,j,k
REAL*8 :: pbar(5),rhobar(5),vxbar(5),vybar(5),vzbar(5),gbar(5),cbar(2)
REAL*8 :: rhojl(2),pjl(2),vxjl(2),vyjl(2),vzjl(2),cjl(2),gjl(2)
REAL*8 :: rhojr(2),pjr(2),vxjr(2),vyjr(2),vzjr(2),cjr(2),gjr(2)
REAL*8 :: deltarho(2),deltap(2),deltavx(2),deltavy(2),deltavz(2),deltac(2),deltag(2)
REAL*8 :: rho6(2),p6(2),vx6(2),vy6(2),vz6(2),c6(2),g6(2)
REAL*8 :: addp,addrho,flat
!
common/bound/pbar,rhobar,vxbar,vybar,vzbar,gbar,cbar
common/jleft/rhojl,pjl,vxjl,vyjl,vzjl,cjl,gjl
common/jright/rhojr,pjr,vxjr,vyjr,vzjr,cjr,gjr
common/delta/deltarho,deltap,deltavx,deltavy,deltavz,deltac,deltag
common/a6/rho6,p6,vx6,vy6,vz6,c6,g6
!
! compute the interpolation functions for the zone 1
! NOTICE that element 1 is not relevant! for cell=1 (=nbound+1=4), elements 2,3,4,5,6 are required

#ifdef STARTING_FROM_1
do j=1,5
   pbar(j)=pres(j+1)
   rhobar(j)=rho(j+1)
   vxbar(j)=vx(j+1)
   vybar(j)=vy(j+1)
   vzbar(j)=vz(j+1)
!   gbar(j)=g(j+1)
enddo
   cbar(1)=c(3)
   cbar(2)=c(4)
#endif
do j=1,5
   pbar(j)=pres(j)
   rhobar(j)=rho(j)
   vxbar(j)=vx(j)
   vybar(j)=vy(j)
   vzbar(j)=vz(j)
!   gbar(j)=g(j)
enddo
   cbar(1)=c(2)
   cbar(2)=c(3)
!
   addp=abs(pbar(4)-pbar(2))/min(pbar(2),pbar(4))
   addrho=0.50*gamma*abs(rhobar(4)-rhobar(2))/&
                  min(rhobar(2),rhobar(4))
   call interp(pbar,pjl,pjr,deltap,p6,flat,0,addp,addrho)
   call interp(rhobar,rhojl,rhojr,deltarho,rho6,flat,1,addp,addrho)
   call interp(vxbar,vxjl,vxjr,deltavx,vx6,flat,0,addp,addrho)
   call interp(gbar,gjl,gjr,deltag,g6,flat,0,addp,addrho)
   call interp(vybar,vyjl,vyjr,deltavy,vy6,flat,0,addp,addrho)
   call interp(vzbar,vzjl,vzjr,deltavz,vz6,flat,0,addp,addrho)
!
END SUBROUTINE boundary

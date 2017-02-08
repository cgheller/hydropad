MODULE ppm_mod
!
!USE dimension
IMPLICIT NONE
SAVE
!
! define basic one dimensional vectors
!
REAL*8  :: rhobar(5),pbar(5),vxbar(5),vybar(5),vzbar(5),gbar(5),cbar(2)
REAL*8  :: rhojl(2),pjl(2),vxjl(2),vyjl(2),vzjl(2),cjl(2),gjl(2)
REAL*8  :: rhojr(2),pjr(2),vxjr(2),vyjr(2),vzjr(2),cjr(2),gjr(2)
REAL*8  :: deltarho(2),deltap(2),deltavx(2),deltavy(2),deltavz(2),deltac(2),deltag(2)
REAL*8  :: rho6(2),p6(2),vx6(2),vy6(2),vz6(2),c6(2),g6(2)
REAL*8  :: pleft,pright,rholeft,rhoright,vxleft,vxright
!$OMP THREADPRIVATE(rhobar,pbar,vxbar,vybar,vzbar,gbar,cbar)
!$OMP THREADPRIVATE(rhojl,pjl,vxjl,vyjl,vzjl,cjl,gjl)
!$OMP THREADPRIVATE(rhojr,pjr,vxjr,vyjr,vzjr,cjr,gjr)
!$OMP THREADPRIVATE(deltarho,deltap,deltavx,deltavy,deltavz,deltac,deltag)
!$OMP THREADPRIVATE(rho6,p6,vx6,vy6,vz6,c6,g6)
!$OMP THREADPRIVATE(pleft,pright,rholeft,rhoright,vxleft,vxright)
!
END MODULE ppm_mod

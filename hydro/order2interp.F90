SUBROUTINE order2interp(ngrid,nmin,nmax,nbound0,pres,rho,vx,vy,vz,g,c,&
           pleft,pright,rholeft,rhoright,vxleft,vxright,vyleft,vyright,vzleft,vzright)
!
!USE dimension
!USE vector
USE scalar
!$acc routine seq
!$acc routine (flatten5) seq
!$acc routine (correction) seq
!$acc routine (beta) seq
!$acc routine (beta0) seq
!
! local variables
!
IMPLICIT NONE
INTEGER, PARAMETER :: n=1
REAL(KIND=8)  :: cbar(2)
REAL(KIND=8)  :: rhojl(2),pjl(2),vxjl(2),vyjl(2),vzjl(2),cjl(2),gjl(2)
REAL(KIND=8)  :: rhojr(2),pjr(2),vxjr(2),vyjr(2),vzjr(2),cjr(2),gjr(2)
REAL(KIND=8)  :: deltarho(2),deltap(2),deltavx(2),deltavy(2),deltavz(2),deltac(2),deltag(2)
REAL(KIND=8)  :: rho6(2),p6(2),vx6(2),vy6(2),vz6(2),c6(2),g6(2)
REAL(KIND=8)  :: pleft,pright,rholeft,rhoright,vxleft,vxright
!
REAL(KIND=8), DIMENSION(ngrid) :: pres,rho,vx,vy
REAL(KIND=8), DIMENSION(ngrid) :: vz,g,c
INTEGER :: mainindex,flatten5
INTEGER :: i,j,k,ngrid,icell,icellaux,cellid,cellidaux,flatyn,ii,icenter
INTEGER :: nmax,nmin,nmid,index,ndir,nbound0
REAL(KIND=8)  :: pcl(3),pcr(3),rhocl(3),rhocr(3),vxcl(3),vxcr(3)
REAL(KIND=8)  :: vycl(3),vycr(3),vzcl(3),vzcr(3),gcl(3),gcr(3)
REAL(KIND=8)  :: ptl,ptr,rhotl,rhotr,vxtl,vxtr,vytl,vytr
REAL(KIND=8)  :: vztl,vztr,ctl,ctr,gtl,gtr
REAL(KIND=8)  :: vyleft,vyright,vzleft,vzright,cleft,cright
REAL(KIND=8)  :: beta,beta0,ee1,ee2,betax
REAL(KIND=8)  :: betaml,beta0l,betapl,betamr,beta0r,betapr
REAL(KIND=8)  :: lambdam,lambdap,lambda0
REAL(KIND=8)  :: flat,addp,addrho,dv2s,dv1s,dvv,dv2r,gradmax,gradleft
!$acc data present(pres,rho,vx,vy,vz,g,c) 
!
! start main loop
!
do i=nmin,nmax
!
! define the starting cell in the stencil at i-th iteration
!
   mainindex = i-nmin+1
!
   flatyn=0
   flat=0.0
   icell = i
   cellidaux = i
#ifdef STENCIL
   cellidaux = cellid
#endif
!CHECK CAREFULLY THIS:
   cbar(1)  = c(i-1)
   cbar(2)  = c(i)
   !cbar(1)  = c(i)
   !cbar(2)  = c(i+1)
!
! calculate flattening parameter
! I had to remove the call to flatten5 and include the loop here because with OpenACC it didn't work
!
   !CALL flatten5(pres(mainindex:mainindex+4),flatyn,5)
   gradmax = 0.0
   do ii = 1,4
      icenter = mainindex+ii
      if(pres(icenter)*pres(icenter-1) .NE. 0.0)&
         gradleft = abs((pres(icenter)-pres(icenter-1))/min(pres(icenter),pres(icenter-1)))
         gradmax = max(gradmax,gradleft)
   enddo
   if(gradmax >= mingradflat)flat=flatvalue
   !if(flatyn == 1)flat=flatvalue
   !CALL flatten5(rho(mainindex:mainindex+4),flatyn,5)
   gradmax = 0.0
   do ii = 1,4
      icenter = mainindex+ii
      if(rho(icenter)*rho(icenter-1) .NE. 0.0)&
         gradleft = abs((rho(icenter)-rho(icenter-1))/min(rho(icenter),rho(icenter-1)))
         gradmax = max(gradmax,gradleft)
   enddo
   if(gradmax >= mingradflat)flat=flatvalue
   !if(flatyn == 1)flat=flatvalue
   !CALL flatten5(vx(mainindex:mainindex+4),flatyn,5)
   gradmax = 0.0
   do ii = 1,4
      icenter = mainindex+ii
      if(vx(icenter)*vx(icenter-1) .NE. 0.0)&
         gradleft = abs((vx(icenter)-vx(icenter-1))/min(vx(icenter),vx(icenter-1)))
         gradmax = max(gradmax,gradleft)
   enddo
   if(gradmax >= mingradflat)flat=flatvalue
   !if(flatyn == 1)flat=flatvalue
   !CALL flatten5(vy(mainindex:mainindex+4),flatyn,5)
   gradmax = 0.0
   do ii = 1,4
      icenter = mainindex+ii
      if(vy(icenter)*vy(icenter-1) .NE. 0.0)&
         gradleft = abs((vy(icenter)-vy(icenter-1))/min(vy(icenter),vy(icenter-1)))
         gradmax = max(gradmax,gradleft)
   enddo
   if(gradmax >= mingradflat)flat=flatvalue
   !if(flatyn == 1)flat=flatvalue
   !CALL flatten5(vz(mainindex:mainindex+4),flatyn,5)
   gradmax = 0.0
   do ii = 1,4
      icenter = mainindex+ii
      if(vz(icenter)*vz(icenter-1) .NE. 0.0)&
         gradleft = abs((vz(icenter)-vz(icenter-1))/min(vz(icenter),vz(icenter-1)))
         gradmax = max(gradmax,gradleft)
   enddo
   if(gradmax >= mingradflat)flat=flatvalue
   !if(flatyn == 1)flat=flatvalue
   !!CALL flattening(flat,ngrid)
   !!flat = 0.0
!
! compute left and right limits of the hydro functions for the
! i-zone
!
   pjl(1)=pjl(2)
   pjr(1)=pjr(2)
   rhojl(1)=rhojl(2)
   rhojr(1)=rhojr(2)
   vxjl(1)=vxjl(2)
   vxjr(1)=vxjr(2)
   vyjl(1)=vyjl(2)
   vyjr(1)=vyjr(2)
   vzjl(1)=vzjl(2)
   vzjr(1)=vzjr(2)
   cjl(1)=cjl(2)
   cjr(1)=cjr(2)
   gjl(1)=gjl(2)
   gjr(1)=gjr(2)
   deltap(1)=deltap(2)
   deltarho(1)=deltarho(2)
   deltavx(1)=deltavx(2)
   deltavy(1)=deltavy(2)
   deltavz(1)=deltavz(2)
   deltac(1)=deltac(2)
   deltag(1)=deltag(2)
   p6(1)=p6(2)
   rho6(1)=rho6(2)
   vx6(1)=vx6(2)
   vy6(1)=vy6(2)
   vz6(1)=vz6(2)
   c6(1)=c6(2)
   g6(1)=g6(2)

! find initial condition for the Riemann problem for each hydro
! variable
!
   addp=abs(pres(i+1)-pres(i-1))/min(pres(i-1),pres(i+1))
   addrho=0.50*gamma*abs(rho(i+1)-rho(i-1))/min(rho(i-1),rho(i+1))
!
   call correction(addp,addrho,pres(mainindex),pjl,pjr,deltap,p6,&
                   cbar,vx(i-1),vx(i),ptl,ptr,pcl,&
                   pcr,flat,0)
   call correction(addp,addrho,rho(mainindex),rhojl,rhojr,deltarho,&
                   rho6,cbar,vx(i-1),vx(i),rhotl,&
                   rhotr,rhocl,rhocr,flat,1)
   call correction(addp,addrho,vx(mainindex),vxjl,vxjr,deltavx,&
                   vx6,cbar,vx(i-1),vx(i),vxtl,&
                   vxtr,vxcl,vxcr,flat,0)
   call correction(addp,addrho,vy(mainindex),vyjl,vyjr,deltavy,&
                   vy6,cbar,vx(i-1),vx(i),vytl,&
                   vytr,vycl,vycr,flat,0)
   call correction(addp,addrho,vz(mainindex),vzjl,vzjr,deltavz,&
                   vz6,cbar,vx(i-1),vx(i),vztl,&
                   vztr,vzcl,vzcr,flat,0)
   call correction(addp,addrho,g(mainindex),gjl,gjr,deltag,&
                   g6,cbar,vx(i-1),vx(i),gtl,&
                   gtr,gcl,gcr,flat,0)
!
           if(ptl.lt.0.0)then
              write(*,*)'negative ptl: ',ptl
              ptl=1e-15
           endif
           if(ptr.lt.0.0)then
              write(*,*)'negative ptr: ',ptr
              ptr=1e-15
           endif
!
           if(rhotl.lt.0.0)then
              write(*,*)'negative rhotl: ',rhotl
              rhotl=1e-15
           endif
           if(rhotr.lt.0.0)then
              write(*,*)'negative rhotr: ',rhotr
              rhotr=1e-15
           endif
!
   ctl=sqrt(gamma*ptl*rhotl)
   ctr=sqrt(gamma*ptr*rhotr)
   lambdam=(vx(i-1)-cbar(1))*rat
   lambda0=vx(i-1)*rat
   lambdap=(vx(i-1)+cbar(1))*rat
   if(lambdam.le.0.0)then
      betaml=0.0
   else
      betax=-1.0
      betaml=beta(ctl,vxtl,vxcl(1),ptl,pcl(1),gcl(1),betax)
   endif           
   if(lambdap.le.0.0)then
      betapl=0.0
      vyleft=vytl
      vzleft=vztl
   else
      betax=1.0
      betapl=beta(ctl,vxtl,vxcl(3),ptl,pcl(3),gcl(3),betax)
      vyleft=vytl
      vzleft=vztl
   endif           
   if(lambda0.le.0.0)then
      beta0l=0.0
   else
      beta0l=beta0(ctl,rhotl,rhocl(2),ptl,pcl(2))
      vyleft=vycl(2)
      vzleft=vzcl(2)
   endif           
!
   lambdam=(vx(i)-cbar(2))*rat
   lambda0=vx(i)*rat
   lambdap=(vx(i)+cbar(2))*rat
   if(lambdap.ge.0.0)then
      betapr=0.0
   else
      betax=1.0
      betapr=beta(ctr,vxtr,vxcr(3),ptr,pcr(3),gcr(3),betax)
   endif           
   if(lambdam.ge.0.0)then
      betamr=0.0
      vyright=vytr
      vzright=vztr
   else
      betax=-1.0
      betamr=beta(ctr,vxtr,vxcr(1),ptr,pcr(1),gcr(1),betax)
      vyright=vytr
      vzright=vztr
   endif           

   if(lambda0.ge.0.0)then
      beta0r=0.0
   else
      beta0r=beta0(ctr,rhotr,rhocr(2),ptr,pcr(2))
      vyright=vycr(2)
      vzright=vzcr(2)
   endif

   pleft=ptl+ctl*ctl*(betapl+betaml)
   if(pleft.le.0.0)pleft=ptl
   pright=ptr+ctr*ctr*(betapr+betamr)
   if(pright.le.0.0)pright=ptr
   vxleft=vxtl+ctl*(betapl-betaml)
   vxright=vxtr+ctr*(betapr-betamr)
   rholeft=1.0/(1.0/rhotl-betaml-beta0l-betapl)
   if(rholeft.le.0.0)rholeft=rhotl
   rhoright=1.0/(1.0/rhotr-betamr-beta0r-betapr)
   if(rhoright.le.0.0)rhoright=rhotr
   cleft=sqrt(gamma*pleft/rholeft)
   cright=sqrt(gamma*pright/rhoright)

enddo   
!$acc end data
!
END SUBROUTINE order2interp

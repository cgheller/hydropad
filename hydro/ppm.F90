SUBROUTINE ppm(ngrid,nbound0,nmin,nmax,cellid,cho,nes,pres,rho,vx,vy,ghalf,vz,g,c,eint,etot,&
               rhoold1d,vxold1d,phi1,phi0,flt,nu,sk2)
!
!USE dimension
!USE vector
USE scalar
!$acc routine seq
!$acc routine (flatten5) seq
!$acc routine (correction) seq
!$acc routine (riemann) seq
!$acc routine (beta) seq
!$acc routine (beta0) seq
!$acc routine (rhoshock) seq
!$acc routine (rarwave) seq
!
! local variables
!
IMPLICIT NONE
INTEGER, PARAMETER :: n=1
REAL*8, DIMENSION(ngrid) :: cho
REAL*8, DIMENSION(ngrid) :: nes
REAL*8  :: rhobar(5),pbar(5),vxbar(5),vybar(5),vzbar(5),gbar(5),cbar(2)
REAL*8  :: rhojl(2),pjl(2),vxjl(2),vyjl(2),vzjl(2),cjl(2),gjl(2)
REAL*8  :: rhojr(2),pjr(2),vxjr(2),vyjr(2),vzjr(2),cjr(2),gjr(2)
REAL*8  :: deltarho(2),deltap(2),deltavx(2),deltavy(2),deltavz(2),deltac(2),deltag(2)
REAL*8  :: rho6(2),p6(2),vx6(2),vy6(2),vz6(2),c6(2),g6(2)
REAL*8  :: pleft,pright,rholeft,rhoright,vxleft,vxright
!
REAL*8, DIMENSION(ngrid) :: pres,rho,vx,vy,ghalf
REAL*8, DIMENSION(ngrid) :: vz,g,c,eint,etot
REAL*8, DIMENSION(ngrid) :: rhoold1d,vxold1d,phi1,phi0
REAL*8, DIMENSION(ngrid) :: flt,nu,sk2
INTEGER :: i,j,k,ngrid,icell,icellaux,cellid,cellidaux,flatyn
INTEGER :: nmax,nmin,nmid,index,ndir,nbound0
REAL*8  :: pnew(ngrid),rhonew(ngrid),vxnew(ngrid)
REAL*8  :: vynew(ngrid),vznew(ngrid),etotnew(ngrid)
REAL*8  :: eintnew(ngrid)
REAL*8  :: flux(2)
REAL*8  :: pcl(3),pcr(3),rhocl(3),rhocr(3),vxcl(3),vxcr(3)
REAL*8  :: vycl(3),vycr(3),vzcl(3),vzcr(3),gcl(3),gcr(3)
REAL*8  :: pm(2),rhom(2),vxm(2),vym(2),vzm(2)
REAL*8  :: ptl,ptr,rhotl,rhotr,vxtl,vxtr,vytl,vytr
REAL*8  :: vztl,vztr,ctl,ctr,gtl,gtr
REAL*8  :: vyleft,vyright,vzleft,vzright,cleft,cright
REAL*8  :: beta,beta0,ee1,ee2,betax
REAL*8  :: betaml,beta0l,betapl,betamr,beta0r,betapr
REAL*8  :: grav,vxold1,lambdam,lambdap,lambda0
REAL*8  :: vvv,dadt
REAL*8  :: flat,addp,addrho,dv2s,dv1s,dvv,dv2r
REAL*8  :: vx1,vx2,p1,p2,rho1,rho2,c1,c2,mmmax
REAL*8  :: rarwave,rhoshock1,rhoshock2,vshock1,vshock2,rhoshock
REAL*8  :: vtail1,vtail2,vhead1,vhead2,rhocont1,rhocont2,vtleft
REAL*8  :: vtright,zeta
REAL*8  :: maxp,pmean
REAL(KIND=8) :: eeleft, ee, eeright, ee0, eeaux
!$acc data present(cho,nes,pres,rho,vx,vy,ghalf,vz,g,c,eint,etot,&
!$acc &            rhoold1d,vxold1d,phi1,phi0,flt,nu,sk2)
!!$acc update device(mingradflat,flatvalue,gamma,rat,m,rgamma1,rm,rgamma,dmax,rdx,rdtath,dt,dx,dat,at,&
!!$acc &              gamma1,gf,eta1,eta2)

!
!CLA
rhoold1d=0.0
vxold1d=0.0
rhonew=0.0
vxnew=0.0
vynew=0.0
vznew=0.0
etotnew=0.0
eintnew=0.0
dadt=1.0-0.5*dat*dt*rat
mmmax=0.0
!#ifndef STENCIL
!flat=flt(3)
!call boundary(flat)
!#endif
!
! nmin and nmax define the actual integration domain: input parameters
! initialize local variables
!
rhom(2)=0.0
vxm(2)=0.0
vym(2)=0.0
vzm(2)=0.0
pm(2)=0.0
!
! start main loop
!
do i=nmin,nmax
   flatyn=0
   flat=0.0
   icell = i
   cellidaux = i
#ifdef STENCIL
   cellidaux = cellid
#endif
!
! select the points involved in the computation for the point i 
!
!CLA
   ghalf=0.0
!CLA

   do j=1,5
      icellaux=icell+j-3
      pbar(j)  = pres(icellaux)
      rhobar(j)= rho(icellaux)
      gbar(j)  = g(icellaux)
      vxbar(j) = vx(icellaux)
      vybar(j) = vy(icellaux)
      vzbar(j) = vz(icellaux)
   enddo
   cbar(1)  = c(icell)
   cbar(2)  = c(icell+1)
   gbar = 0.0
!
! calculate flattening parameter
!
   CALL flatten5(pbar,flatyn,5)
   if(flatyn == 1)flat=flatvalue
   CALL flatten5(rhobar,flatyn,5)
   if(flatyn == 1)flat=flatvalue
   CALL flatten5(vxbar,flatyn,5)
   if(flatyn == 1)flat=flatvalue
   CALL flatten5(vybar,flatyn,5)
   if(flatyn == 1)flat=flatvalue
   CALL flatten5(vzbar,flatyn,5)
   if(flatyn == 1)flat=flatvalue
   
#ifdef PIPPO
!PIPPO
#endif
   !CALL flattening(flat,ngrid)
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
   addp=abs(pbar(4)-pbar(2))/min(pbar(2),pbar(4))
   addrho=0.50*gamma*abs(rhobar(4)-rhobar(2))/&
                           min(rhobar(2),rhobar(4))
   call correction(addp,addrho,pbar,pjl,pjr,deltap,p6,&
                   cbar,vxbar(2),vxbar(3),ptl,ptr,pcl,&
                   pcr,flat,0)
   call correction(addp,addrho,rhobar,rhojl,rhojr,deltarho,&
                   rho6,cbar,vxbar(2),vxbar(3),rhotl,&
                   rhotr,rhocl,rhocr,flat,1)
   call correction(addp,addrho,vxbar,vxjl,vxjr,deltavx,&
                   vx6,cbar,vxbar(2),vxbar(3),vxtl,&
                   vxtr,vxcl,vxcr,flat,0)
   call correction(addp,addrho,vybar,vyjl,vyjr,deltavy,&
                   vy6,cbar,vxbar(2),vxbar(3),vytl,&
                   vytr,vycl,vycr,flat,0)
   call correction(addp,addrho,vzbar,vzjl,vzjr,deltavz,&
                   vz6,cbar,vxbar(2),vxbar(3),vztl,&
                   vztr,vzcl,vzcr,flat,0)
   call correction(addp,addrho,gbar,gjl,gjr,deltag,&
                   g6,cbar,vxbar(2),vxbar(3),gtl,&
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
   lambdam=(vxbar(2)-cbar(1))*rat
   lambda0=vxbar(2)*rat
   lambdap=(vxbar(2)+cbar(1))*rat
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
   lambdam=(vxbar(3)-cbar(2))*rat
   lambda0=vxbar(3)*rat
   lambdap=(vxbar(3)+cbar(2))*rat
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

   IF (i.eq.nmin) CYCLE

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
   if(pleft.le.pright)then
     p1=pleft
     p2=pright
     vx1=vxleft
     vx2=vxright
     rho1=rholeft
     rho2=rhoright
     c1=cleft
     c2=cright
   else
     p1=pright
     p2=pleft
     vx1=vxleft
     vx2=vxright
     rho1=rhoright
     rho2=rholeft
     c1=cright
     c2=cleft
   endif
   dvv=vx1-vx2
   dv2s=(p2-p1)/sqrt(0.50*rho1*((gamma-1)*p1+(gamma+1)*p2))
   dv1s=-(2.0*c2/(gamma-1.0))*(1.0-(p1/p2)**m)
   dv2r=-2.0*(c1+c2)*rgamma1
!
! solve the Riemann problem ( 1=(i-1/2), 2=(i+1/2) ).
! (the suffix m indicate the time averaged quantities at the zone boundary)
! for (i-1/2):
!
   pm(1)=pm(2)
   rhom(1)=rhom(2)
   vxm(1)=vxm(2)
   vym(1)=vym(2)
   vzm(1)=vzm(2)
   !write(*,*)i,rhom(1),rhom(2)
!
! for (i+1/2) :
! 1 > two shocks
!
   if(dvv.ge.dv2s)then
!
! compute post-shock pressure and velocity (same for both the shocks):
!
      pm(2)=max(pleft,pright)
      call riemann(pm(2),vxm(2),pleft,pright,rholeft,rhoright,vxleft,vxright)
!
! compute post-shock density and velocity for left (1) and right (2) shock:
!
      rhoshock1=rhoshock(pm(2),pleft,rholeft)
      rhoshock2=rhoshock(pm(2),pright,rhoright)
      if(rhoshock1.eq.rholeft)then
            vshock1=0.0
      else
            vshock1=(rhoshock1*vxm(2)-rholeft*vxleft)/(rhoshock1-rholeft)
      endif
      if(rhoshock2.eq.rhoright)then
            vshock2=0.0
      else
            vshock2=(rhoshock2*vxm(2)-rhoright*vxright)/(rhoshock2-rhoright)
      endif
!
! solve Riemann problem:
!
      if(vshock2.lt.0.0)then
         pm(2)=pright
         vxm(2)=vxright
         rhom(2)=rhoright
      else if(vxm(2).le.0.0.and.vshock2.ge.0.0)then
         rhom(2)=rhoshock2
      else if(vshock1.le.0.0.and.vxm(2).gt.0.0)then
         rhom(2)=rhoshock1
      else
         pm(2)=pleft
         vxm(2)=vxleft
         rhom(2)=rholeft
      endif
!
! 2 > one shock and one rarefaction wave
!
   else if(dvv.lt.dv2s.and.dvv.ge.dv1s)then
!
! compute post-shock pressure and velocity:
!
      pm(2)=( (cleft+cright+(vxleft-vxright)*&
                      (gamma-1)*0.5)/(cleft/pleft**m+cright/&
                       pright**m) )**rm
              call riemann(pm(2),vxm(2),pleft,pright,rholeft,rhoright,vxleft,vxright)
!
! solve Riemann problem in the case pleft > pright:
!
      if(pleft.ge.pright)then
	 rhoshock2=rhoshock(pm(2),pright,rhoright)
	 rhocont1=rholeft*(pm(2)/pleft)**rgamma
	 if(rhoshock2.eq.rhoright)then
	    vshock2=0.0
	 else

	    vshock2=(rhoshock2*vxm(2)-rhoright*vxright)/&
                            (rhoshock2-rhoright)
	 endif
	 vtail1=vxm(2)-sqrt(gamma*pm(2)/rhocont1)
	 vhead1=vxleft-cleft
	 if(vshock2.lt.0)then
	    pm(2)=pright
	    vxm(2)=vxright
	    rhom(2)=rhoright
	 else if(vhead1.ge.0.0)then
	    pm(2)=pleft
	    vxm(2)=vxleft
	    rhom(2)=rholeft
	 else if(vxm(2).le.0.0.and.vshock2.ge.0)then
	    rhom(2)=rhoshock2
	 else if(vxm(2).gt.0.0.and.vtail1.le.0.0)then
	    rhom(2)=rhocont1
	 else if(vtail1.gt.0.0.and.vhead1.le.0.0)then
	    pm(2)=rarwave(pleft,pm(2),vtail1,vhead1)
	    rhom(2)=rarwave(rholeft,rhocont1,vtail1,vhead1)
	    vxm(2)=rarwave(vxleft,vxm(2),vtail1,vhead1)
	 endif
!
! solve Riemann problem in the case pleft < pright:
!
      else       
	 rhoshock1=rhoshock(pm(2),pleft,rholeft)
	 rhocont2=rhoright*(pm(2)/pright)**rgamma
	 if(rhoshock1.eq.rholeft)then
	    vshock1=0.0
	 else
	    vshock1=(rhoshock1*vxm(2)-rholeft*vxleft)/&
                            (rhoshock1-rholeft)
	 endif
	 vtail2=vxm(2)+sqrt(gamma*pm(2)/rhocont2)
	 vhead2=vxright+cright
	 if(vshock1.gt.0)then
	    pm(2)=pleft
	    vxm(2)=vxleft
	    rhom(2)=rholeft
	 else if(vhead2.le.0.0)then
	    pm(2)=pright
	    vxm(2)=vxright
	    rhom(2)=rhoright
	 else if(vxm(2).ge.0.0.and.vshock1.le.0.0)then
	    rhom(2)=rhoshock1
	 else if(vxm(2).lt.0.0.and.vtail2.ge.0.0)then
	    rhom(2)=rhocont2
	 else if(vtail2.lt.0.0.and.vhead2.ge.0.0)then
	    pm(2)=rarwave(pright,pm(2),vtail2,vhead2)
	    rhom(2)=rarwave(rhoright,rhocont2,vtail2,vhead2)
	    vxm(2)=rarwave(vxright,vxm(2),vtail2,vhead2)
	 endif
      endif
!
! 3 > two rarefaction wave
!
   else if(dvv.lt.dv1s.and.dvv.gt.dv2r)then
!
! compute pressure and velocity behind the rarefactions wave:
!
      zeta=(cright/cleft)*(pleft/pright)**m
      vtleft=vxleft+2.0*cleft*rgamma1
      vtright=vxright-2.0*cright*rgamma1
      vxm(2)=(vtleft*zeta+vtright)/(1.0+zeta)
      pm(2)=pleft*((gamma-1)*(vtleft-vtright)/&
                    (2.0*cleft*(1.0+zeta)))**rm
      rhocont1=rholeft*(pm(2)/pleft)**rgamma
      rhocont2=rhoright*(pm(2)/pright)**rgamma
      vtail2=vxm(2)+sqrt(gamma*pm(2)/rhocont2)
      vhead2=vxright+cright
      vtail1=vxm(2)-sqrt(gamma*pm(2)/rhocont1)
      vhead1=vxleft-cleft
!
! solve Riemann problem:
!
      if(vhead2.le.0.0)then
	    pm(2)=pright
	    vxm(2)=vxright
	    rhom(2)=rhoright
      else if(vhead1.ge.0.0)then
	    pm(2)=pleft
	    vxm(2)=vxleft
	    rhom(2)=rholeft
      else if(vhead2.gt.0.0.and.vtail2.le.0.0)then
	    pm(2)=rarwave(pright,pm(2),vtail2,vhead2)
	    rhom(2)=rarwave(rhoright,rhocont2,vtail2,vhead2)
	    vxm(2)=rarwave(vxright,vxm(2),vtail2,vhead2)
      else if(vxm(2).le.0.0.and.vtail2.gt.0.0)then
	    rhom(2)=rhocont2
      else if(vxm(2).gt.0.0.and.vtail1.le.0.0)then
	    rhom(2)=rhocont1
      else if(vtail1.gt.0.0.and.vhead1.le.0.0)then
	    pm(2)=rarwave(pleft,pm(2),vtail1,vhead1)
	    rhom(2)=rarwave(rholeft,rhocont1,vtail1,vhead1)
	    vxm(2)=rarwave(vxleft,vxm(2),vtail1,vhead1)
      endif
!
! 4 > vacuum zone
!
   else if(dvv.le.dv2r)then
!
! compute pressure and velocity behind the rarefactions waves:
!
      pm(2)=0.0
      vxm(2)=0.0
      rhocont1=0.0
      rhocont2=0.0
      vtail2=vxright-2.0*cright*rgamma1
      vhead2=vxright+cright
      vtail1=vxleft-2.0*cleft*rgamma1
      vhead1=vxleft-cleft
!
! solve Riemann problem:
!
      if(vhead2.le.0.0)then
	    pm(2)=pright
	    vxm(2)=vxright
	    rhom(2)=rhoright
      else if(vhead1.ge.0.0)then
	    pm(2)=pleft
	    vxm(2)=vxleft
	    rhom(2)=rholeft
      else if(vhead2.gt.0.0.and.vtail2.le.0.0)then
	    pm(2)=rarwave(pright,pm(2),vtail2,vhead2)
	    rhom(2)=rarwave(rhoright,rhocont2,vtail2,vhead2)
	    vxm(2)=rarwave(vxright,vxm(2),vtail2,vhead2)
      else if(vtail1.le.0.0.and.vtail2.gt.0.0)then
	    rhom(2)=0.0
      else if(vtail1.gt.0.0.and.vhead1.le.0.0)then
	    pm(2)=rarwave(pleft,pm(2),vtail1,vhead1)
	    rhom(2)=rarwave(rholeft,rhocont1,vtail1,vhead1)
	    vxm(2)=rarwave(vxleft,vxm(2),vtail1,vhead1)
      endif
   endif
!
! compute j+1/2 value for y and z velocity components
!
   if(vxm(2).ge.0.0)then
      vym(2)=dadt*vyleft
      vzm(2)=dadt*vzleft
   else
      vym(2)=dadt*vyright
      vzm(2)=dadt*vzright
   endif
!
! correct pressure on very high density peaks
!
   if(rhobar(2).gt.dmax*rhobar(1).and.&
              rhobar(3).gt.dmax*rhobar(4))then
                  maxp=(vxbar(2)-vxbar(3))*&
                       max(rhobar(2),rhobar(3))*&
                       max(mmmax,vxbar(2)-vxbar(3))
                  pm(2)=max(pm(2),maxp)
   endif
!
!CLA

   nmid=nmin+1
   IF (i.eq.nmid) CYCLE
   !write(*,*)icell
   !write(*,*)pm(1),rhom(1),vxm(1)
   !write(*,*)pm(2),rhom(2),vxm(2)
!
! solve the hydro equations for rho, vx, vy, vz, p (note that
! for the variable flux 1=i-(1/2), 2=i+(1/2)):
! rho >
!
   index=icell-1
!
   flux(1)=rhom(1)*vxm(1)
   flux(2)=rhom(2)*vxm(2)
   rhonew(index)=rho(index)+rdtath*(flux(1)-flux(2))*rdx
   !write(*,*)index,icell
   !write(*,*)flux(1),flux(2),rhonew(index)
   !write(*,*)"--------------------------"
! vx >
   vxold1=vx(index)
   flux(1)=rhom(1)*vxm(1)*vxm(1)+pm(1)
   flux(2)=rhom(2)*vxm(2)*vxm(2)+pm(2)
   grav=(rhoold1d(index)+rhonew(index))*ghalf(index)*0.5
   vxnew(index)=(vx(index)*rho(index)+&
                         rdtath*(flux(1)-flux(2))*rdx+rdtath*grav)&
                         /rhonew(index)
! vy >
   flux(1)=vym(1)*vxm(1)*rhom(1)
   flux(2)=vym(2)*vxm(2)*rhom(2)
   vynew(index)=(rho(index)*&
                        vy(index)+rdtath*(flux(1)-flux(2))*rdx)&
                        /rhonew(index)
! vz >
   flux(1)=vzm(1)*vxm(1)*rhom(1)
   flux(2)=vzm(2)*vxm(2)*rhom(2)
   vznew(index)=(rho(index)*&
                        vz(index)+rdtath*(flux(1)-flux(2))*rdx)&
                        /rhonew(index)
! Etot >
   flux(1)=(gamma*pm(1)*rgamma1+0.50*rhom(1)*(vxm(1)*vxm(1)+vym(1)*vym(1)+vzm(1)*vzm(1)))*vxm(1)
   flux(2)=(gamma*pm(2)*rgamma1+0.50*rhom(2)*(vxm(2)*vxm(2)+vym(2)*vym(2)+vzm(2)*vzm(2)))*vxm(2)
   grav=(rhoold1d(index)*vxold1d(index)+&
                 rhonew(index)*vxnew(index))&
                 *ghalf(index)*0.5
   etotnew(index)=etot(index)+rdtath*(flux(1)-flux(2))*rdx+rdtath*grav
! Eint >
! CLA
! THE EXPANSION TERM SEEMS TO BE MISSING
   pmean=(pm(1)+pm(2))*0.5
   if(rhobar(2).gt.dmax*rhobar(1).and.&
      rhobar(2).gt.dmax*rhobar(3))then
            if(cho(index).eq.1.0)then
               maxp=(vxm(1)-vxm(2))*rhobar(2)*max(mmmax,vxm(1)-vxm(2))
               pmean=max(pmean,maxp)
            else if(cho(index).eq.2.0)then
               pmean=(gamma-1.0)*eint(index)
            endif
   endif
   eintnew(index)=eint(index)+rdtath*pmean*(vxm(1)-vxm(2))*rdx
   flux(1)=(pm(1)*vxm(1))*rgamma1
   flux(2)=(pm(2)*vxm(2))*rgamma1
   eintnew(index)=eintnew(index)+rdtath*(flux(1)-flux(2))*rdx
   if(eintnew(index).lt.0)then
         eintnew(index)=eint(index)
   endif
!
! end of the main loop
!
enddo
!
! use only internal energy or both total and internal energy 
!
#ifndef STENCIL
   do i=nmin,nmax
#else
   i=4
#endif
!CLA this introduces differences in parallel results
      eeleft = etotnew(i-1)
      if(i .EQ. nmin) eeleft = 0.0
      ee = etotnew(i)
      eeright = etotnew(i+1)
      if(i .EQ. nmax) eeright = 0.0
      ee2 = max(eeleft, ee, eeright)
#ifdef STENCIL
      ee2 = ee
#endif
      vvv = vxnew(i)*vxnew(i)+vynew(i)*vynew(i)+vznew(i)*vznew(i)
      eeaux = etotnew(i)-0.50*rhonew(i)*vvv 
      ee1 = eeaux/ee2
      if (ee1 >= eta2) eintnew(i)=eeaux
       
      ee0 = eeaux/etotnew(i)
      if (ee0 >= eta1) then
         pnew(i) = (gamma-1)*eeaux
      else
         pnew(i) = (gamma-1)*eintnew(i)
         !write(*,*)"Correcting with internal energy!!!!!!!!!!!!!!!"
      endif   
#ifndef STENCIL
   enddo
#endif
!
! update of the main variables
!
#ifndef STENCIL
   do i=nmin,nmax
#else
   i=4
#endif
  if(eintnew(i).le.0.0)then
     eintnew(i)=eint(i)
     pnew(i)=pres(i)
  endif
  if(rhonew(i).le.0.)rhonew(i)=rho(i)
  pres(i)=pnew(i)  
  rho(i)=rhonew(i)      
  vx(i)=vxnew(i)        
  vy(i)=vynew(i)        
  vz(i)=vznew(i)        
  etot(i)=etotnew(i)    
  eint(i)=eintnew(i)
#ifndef STENCIL
enddo   
#endif
!$acc end data
!
END SUBROUTINE ppm

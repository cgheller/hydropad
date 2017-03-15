SUBROUTINE riemannsolver(ngrid,rhobar,vxbar,pleft,pright,rholeft,rhoright,vxleft,vxright,vyleft,&
                         vyright,vzleft,vzright,pm,rhom,vxm,vym,vzm)
!
!USE dimension
!USE vector
USE scalar
!$acc routine seq
!$acc routine (rhoshock) seq
!$acc routine (riemann) seq
!$acc routine (rarwave) seq
!
! local variables
!
IMPLICIT NONE
INTEGER, PARAMETER :: n=1
REAL(KIND=8)  :: pleft,pright,rholeft,rhoright,vxleft,vxright
REAL(KIND=8)  :: vyleft,vyright,vzleft,vzright,cleft,cright
!
INTEGER :: i,j,k,ngrid,icell,icellaux,cellid,cellidaux,flatyn
INTEGER :: nmax,nmin,nmid,index,ndir,nbound0
REAL(KIND=8)  :: rhobar(ngrid),vxbar(ngrid)
REAL(KIND=8)  :: pm,rhom,vxm,vym,vzm
REAL(KIND=8)  :: dadt
REAL(KIND=8)  :: flat,addp,addrho,dv2s,dv1s,dvv,dv2r
REAL(KIND=8)  :: vx1,vx2,p1,p2,rho1,rho2,c1,c2,mmmax
REAL(KIND=8)  :: rarwave,rhoshock1,rhoshock2,vshock1,vshock2,rhoshock
REAL(KIND=8)  :: vtail1,vtail2,vhead1,vhead2,rhocont1,rhocont2,vtleft
REAL(KIND=8)  :: vtright,zeta
REAL(KIND=8)  :: maxp,pmean
!$acc data create(rhobar,vxbar)

cleft  = sqrt(gamma*pleft/rholeft)
cright = sqrt(gamma*pright/rhoright)
!
!CLA
dadt=1.0-0.5*dat*dt*rat
mmmax=0.0
!
! solve the Riemann problem
!
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
!
dvv=vx1-vx2
dv2s=(p2-p1)/sqrt(0.50*rho1*((gamma-1)*p1+(gamma+1)*p2))
dv1s=-(2.0*c2/(gamma-1.0))*(1.0-(p1/p2)**m)
dv2r=-2.0*(c1+c2)*rgamma1
!
! solve the Riemann problem ( 1=(i-1/2), 2=(i+1/2) ).
! (the suffix m indicate the time averaged quantities at the zone boundary)
! for (i-1/2):
!
! for (i-1/2) :
! 1 > two shocks
!
if(dvv.ge.dv2s)then
!
! compute post-shock pressure and velocity (same for both the shocks):
!
   pm=max(pleft,pright)
   call riemann(pm,vxm,pleft,pright,rholeft,rhoright,vxleft,vxright)
!
! compute post-shock density and velocity for left (1) and right (2) shock:
!
   rhoshock1=rhoshock(pm,pleft,rholeft)
   rhoshock2=rhoshock(pm,pright,rhoright)
   if(rhoshock1.eq.rholeft)then
            vshock1=0.0
   else
            vshock1=(rhoshock1*vxm-rholeft*vxleft)/(rhoshock1-rholeft)
   endif
   if(rhoshock2.eq.rhoright)then
            vshock2=0.0
   else
            vshock2=(rhoshock2*vxm-rhoright*vxright)/(rhoshock2-rhoright)
   endif
!
! solve Riemann problem:
!
   if(vshock2.lt.0.0)then
         pm=pright
         vxm=vxright
         rhom=rhoright
   else if(vxm.le.0.0.and.vshock2.ge.0.0)then
         rhom=rhoshock2
   else if(vshock1.le.0.0.and.vxm.gt.0.0)then
         rhom=rhoshock1
   else
         pm=pleft
         vxm=vxleft
         rhom=rholeft
   endif
!
! 2 > one shock and one rarefaction wave
!
else if(dvv.lt.dv2s.and.dvv.ge.dv1s)then
!
! compute post-shock pressure and velocity:
!
  pm=( (cleft+cright+(vxleft-vxright)*&
                   (gamma-1)*0.5)/(cleft/pleft**m+cright/&
                   pright**m) )**rm
  call riemann(pm,vxm,pleft,pright,rholeft,rhoright,vxleft,vxright)
!
! solve Riemann problem in the case pleft > pright:
!
  if(pleft.ge.pright)then
         rhoshock2=rhoshock(pm,pright,rhoright)
         rhocont1=rholeft*(pm/pleft)**rgamma
         if(rhoshock2.eq.rhoright)then
          vshock2=0.0
         else

         vshock2=(rhoshock2*vxm-rhoright*vxright)/&
                            (rhoshock2-rhoright)
      endif
      vtail1=vxm-sqrt(gamma*pm/rhocont1)
      vhead1=vxleft-cleft
      if(vshock2.lt.0)then
        pm=pright
        vxm=vxright
        rhom=rhoright
      else if(vhead1.ge.0.0)then
        pm=pleft
        vxm=vxleft
        rhom=rholeft
      else if(vxm.le.0.0.and.vshock2.ge.0)then
        rhom=rhoshock2
      else if(vxm.gt.0.0.and.vtail1.le.0.0)then
        rhom=rhocont1
      else if(vtail1.gt.0.0.and.vhead1.le.0.0)then
        pm=rarwave(pleft,pm,vtail1,vhead1)
        rhom=rarwave(rholeft,rhocont1,vtail1,vhead1)
        vxm=rarwave(vxleft,vxm,vtail1,vhead1)
      endif
!
! solve Riemann problem in the case pleft < pright:
!
  else       
      rhoshock1=rhoshock(pm,pleft,rholeft)
      rhocont2=rhoright*(pm/pright)**rgamma
      if(rhoshock1.eq.rholeft)then
        vshock1=0.0
      else
        vshock1=(rhoshock1*vxm-rholeft*vxleft)/&
                              (rhoshock1-rholeft)
      endif
        vtail2=vxm+sqrt(gamma*pm/rhocont2)
        vhead2=vxright+cright
      if(vshock1.gt.0)then
        pm=pleft
        vxm=vxleft
        rhom=rholeft
      else if(vhead2.le.0.0)then
        pm=pright
        vxm=vxright
        rhom=rhoright
      else if(vxm.ge.0.0.and.vshock1.le.0.0)then
        rhom=rhoshock1
      else if(vxm.lt.0.0.and.vtail2.ge.0.0)then
        rhom=rhocont2
      else if(vtail2.lt.0.0.and.vhead2.ge.0.0)then
        pm=rarwave(pright,pm,vtail2,vhead2)
        rhom=rarwave(rhoright,rhocont2,vtail2,vhead2)
        vxm=rarwave(vxright,vxm,vtail2,vhead2)
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
  vxm=(vtleft*zeta+vtright)/(1.0+zeta)
  pm=pleft*((gamma-1)*(vtleft-vtright)/&
                (2.0*cleft*(1.0+zeta)))**rm
  rhocont1=rholeft*(pm/pleft)**rgamma
  rhocont2=rhoright*(pm/pright)**rgamma
  vtail2=vxm+sqrt(gamma*pm/rhocont2)
  vhead2=vxright+cright
  vtail1=vxm-sqrt(gamma*pm/rhocont1)
  vhead1=vxleft-cleft
!
! solve Riemann problem:
!
  if(vhead2.le.0.0)then
    pm=pright
    vxm=vxright
    rhom=rhoright
  else if(vhead1.ge.0.0)then
    pm=pleft
    vxm=vxleft
    rhom=rholeft
  else if(vhead2.gt.0.0.and.vtail2.le.0.0)then
    pm=rarwave(pright,pm,vtail2,vhead2)
    rhom=rarwave(rhoright,rhocont2,vtail2,vhead2)
    vxm=rarwave(vxright,vxm,vtail2,vhead2)
  else if(vxm.le.0.0.and.vtail2.gt.0.0)then
    rhom=rhocont2
      else if(vxm.gt.0.0.and.vtail1.le.0.0)then
    rhom=rhocont1
  else if(vtail1.gt.0.0.and.vhead1.le.0.0)then
    pm=rarwave(pleft,pm,vtail1,vhead1)
    rhom=rarwave(rholeft,rhocont1,vtail1,vhead1)
    vxm=rarwave(vxleft,vxm,vtail1,vhead1)
  endif
!
! 4 > vacuum zone
!
else if(dvv.le.dv2r)then
!
! compute pressure and velocity behind the rarefactions waves:
!
   pm=0.0
   vxm=0.0
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
    pm=pright
    vxm=vxright
    rhom=rhoright
  else if(vhead1.ge.0.0)then
    pm=pleft
    vxm=vxleft
    rhom=rholeft
  else if(vhead2.gt.0.0.and.vtail2.le.0.0)then
    pm=rarwave(pright,pm,vtail2,vhead2)
    rhom=rarwave(rhoright,rhocont2,vtail2,vhead2)
    vxm=rarwave(vxright,vxm,vtail2,vhead2)
  else if(vtail1.le.0.0.and.vtail2.gt.0.0)then
    rhom=0.0
  else if(vtail1.gt.0.0.and.vhead1.le.0.0)then
    pm=rarwave(pleft,pm,vtail1,vhead1)
    rhom=rarwave(rholeft,rhocont1,vtail1,vhead1)
    vxm=rarwave(vxleft,vxm,vtail1,vhead1)
  endif
endif
!
! compute j+1/2 value for y and z velocity components
!
if(vxm.ge.0.0)then
      vym=dadt*vyleft
      vzm=dadt*vzleft
else
      vym=dadt*vyright
      vzm=dadt*vzright
endif
!
! correct pressure on very high density peaks
!
if(rhobar(2).gt.dmax*rhobar(1).and.&
   rhobar(3).gt.dmax*rhobar(4))then
              maxp=(vxbar(2)-vxbar(3))*&
                   max(rhobar(2),rhobar(3))*&
                   max(mmmax,vxbar(2)-vxbar(3))
              pm=max(pm,maxp)
endif
!
!$acc end data
!
END SUBROUTINE riemannsolver

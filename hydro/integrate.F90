SUBROUTINE integrate(pm1,pm2,rhom1,rhom2,vxm1,vxm2,vym1,vym2,vzm1,vzm2,&
                     pres,rho,vx,vy,vz,etot,eint,ghalf)
!
!USE dimension
!USE vector
USE scalar
!$acc routine seq
!
! local variables
!
IMPLICIT NONE
REAL(KIND=8)  :: pnew,rhonew,vxnew,vynew,vznew,etotnew,eintnew,ghalf
REAL(KIND=8)  :: pres,rho,vx,vy,vz,etot,eint
REAL(KIND=8)  :: vxold1d
REAL(KIND=8)  :: flux1,flux2
REAL(KIND=8)  :: pm1,rhom1,vxm1,vym1,vzm1
REAL(KIND=8)  :: pm2,rhom2,vxm2,vym2,vzm2
REAL(KIND=8)  :: vvv,dadt
REAL(KIND=8)  :: mmmax,grav,pmean,ee1,ee2
REAL(KIND=8) :: eeleft, ee, eeright, ee0, eeaux

vxold1d=0.0
rhonew=0.0
vxnew=0.0
vynew=0.0
vznew=0.0
etotnew=0.0
eintnew=0.0
dadt=1.0-0.5*dat*dt*rat
mmmax=0.0
!
! solve the hydro equations for rho, vx, vy, vz, p (note that
! for the variable flux 1=i-(1/2), 2=i+(1/2)):
! rho >
!
   flux1=rhom1*vxm1
   flux2=rhom2*vxm2
   rhonew=rho+rdtath*(flux1-flux2)*rdx
! vx >
   vxold1d=vx
   flux1=rhom1*vxm1*vxm1+pm1
   flux2=rhom2*vxm2*vxm2+pm2
   grav=(rho+rhonew)*ghalf*0.5
   vxnew=(vx*rho+rdtath*(flux1-flux2)*rdx+rdtath*grav)/rhonew
! vy >
   flux1=vym1*vxm1*rhom1
   flux2=vym2*vxm2*rhom2
   vynew=(rho*vy+rdtath*(flux1-flux2)*rdx)/rhonew
! vz >
   flux1=vzm1*vxm1*rhom1
   flux2=vzm2*vxm2*rhom2
   vznew=(rho*vz+rdtath*(flux1-flux2)*rdx)/rhonew
! Etot >
   flux1=(gamma*pm1*rgamma1+0.50*rhom1*(vxm1*vxm1+vym1*vym1+vzm1*vzm1))*vxm1
   flux2=(gamma*pm2*rgamma1+0.50*rhom2*(vxm2*vxm2+vym2*vym2+vzm2*vzm2))*vxm2
   grav=(rho*vx+rhonew*vxnew)*ghalf*0.5
   etotnew=etot+rdtath*(flux1-flux2)*rdx+rdtath*grav
! Eint >
! I don't understand this:
#ifdef HIGHPEAKS
   pmean=(pm1+pm2)*0.5
   if(rhobar(2).gt.dmax*rhobar(1).and.&
      rhobar(2).gt.dmax*rhobar(3))then
            if(cho.eq.1.0)then
               maxp=(vxm1-vxm2)*rhobar(2)*max(mmmax,vxm1-vxm2)
               pmean=max(pmean,maxp)
            else if(cho.eq.2.0)then
               pmean=(gamma-1.0)*eint
            endif
   endif
#endif
!
! CLA
! THE EXPANSION TERM SEEMS TO BE MISSING
   eintnew=eint+rdtath*pmean*(vxm1-vxm2)*rdx
   flux1=(pm1*vxm1)*rgamma1
   flux2=(pm2*vxm2)*rgamma1
   eintnew=eintnew+rdtath*(flux1-flux2)*rdx
   if(eintnew.lt.0)then
         eintnew=eint
   endif
! use only internal energy or both total and internal energy 
!
!CLA this introduces differences in parallel results
  ee = etotnew
  ee2 = ee
  vvv = vxnew*vxnew+vynew*vynew+vznew*vznew
  eeaux = etotnew-0.50*rhonew*vvv 
  ee1 = eeaux/ee2
  if (ee1 >= eta2) eintnew=eeaux
 
  ee0 = eeaux/etotnew
  if (ee0 >= eta1) then
     pnew = (gamma-1)*eeaux
  else
     pnew = (gamma-1)*eintnew
     !write(*,*)"Correcting with internal energy!!!!!!!!!!!!!!!"
  endif   
!
! update of the main variables
!
  if(eintnew.le.0.0)then
     eintnew=eint
     pnew=pres
  endif
  if(rhonew.le.0.)rhonew=rho
  pres=pnew  
  rho=rhonew      
  vx=vxnew        
  vy=vynew        
  vz=vznew        
!
END SUBROUTINE integrate

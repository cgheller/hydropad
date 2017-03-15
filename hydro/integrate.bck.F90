SUBROUTINE integrate(pm(1),pm(2),rhom(1),rhom(2),vxm(1),vxm(2),vym(1),vym(2),vzm(1),vzm(2),&
                     pres,rho,vx,vy,vz,etot,eint)
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
REAL(KIND=8)  :: vxold1d,rhoold1d
REAL(KIND=8)  :: flux(2)
REAL(KIND=8)  :: pm(2),rhom(2),vxm(2),vym(2),vzm(2)
REAL(KIND=8)  :: vvv,dadt
REAL(KIND=8) :: eeleft, ee, eeright, ee0, eeaux

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
ghalf=0.0
!
! solve the hydro equations for rho, vx, vy, vz, p (note that
! for the variable flux 1=i-(1/2), 2=i+(1/2)):
! rho >
!
   flux(1)=rhom(1)*vxm(1)
   flux(2)=rhom(2)*vxm(2)
   rhonew=rho+rdtath*(flux(1)-flux(2))*rdx
! vx >
   vxold1d=vx
   flux(1)=rhom(1)*vxm(1)*vxm(1)+pm(1)
   flux(2)=rhom(2)*vxm(2)*vxm(2)+pm(2)
   grav=(rho+rhonew)*ghalf*0.5
   vxnew=(vx*rho+rdtath*(flux(1)-flux(2))*rdx+rdtath*grav)/rhonew(index)
! vy >
   flux(1)=vym(1)*vxm(1)*rhom(1)
   flux(2)=vym(2)*vxm(2)*rhom(2)
   vynew=(rho*vy+rdtath*(flux(1)-flux(2))*rdx)/rhonew
! vz >
   flux(1)=vzm(1)*vxm(1)*rhom(1)
   flux(2)=vzm(2)*vxm(2)*rhom(2)
   vznew=(rho*vz+rdtath*(flux(1)-flux(2))*rdx)/rhonew(index)
! Etot >
   flux(1)=(gamma*pm(1)*rgamma1+0.50*rhom(1)*(vxm(1)*vxm(1)+vym(1)*vym(1)+vzm(1)*vzm(1)))*vxm(1)
   flux(2)=(gamma*pm(2)*rgamma1+0.50*rhom(2)*(vxm(2)*vxm(2)+vym(2)*vym(2)+vzm(2)*vzm(2)))*vxm(2)
   grav=(rho*vx+rhonew*vxnew)*ghalf*0.5
   etotnew=etot+rdtath*(flux(1)-flux(2))*rdx+rdtath*grav
! Eint >
! CLA
! THE EXPANSION TERM SEEMS TO BE MISSING
#ifdef HIGHPEAKS
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
#endif
   eintnew=eint+rdtath*pmean*(vxm(1)-vxm(2))*rdx
   flux(1)=(pm(1)*vxm(1))*rgamma1
   flux(2)=(pm(2)*vxm(2))*rgamma1
   eintnew=eintnew+rdtath*(flux(1)-flux(2))*rdx
   if(eintnew.lt.0)then
         eintnew=eint
   endif
!
! end of the main loop
!
enddo
!
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
  if(rhonew(i).le.0.)rhonew(i)=rho(i)
  pres=pnew  
  rho=rhonew      
  vx=vxnew        
  vy=vynew        
  vz=vznew        
#ifndef STENCIL
enddo   
#endif
!$acc end data
!
END SUBROUTINE integrate

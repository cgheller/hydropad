SUBROUTINE flattening(flat,ngrid,j)
!
USE vector
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k,indexj
INTEGER :: ngrid,nbound
INTEGER :: sj,sj1
REAL*8  :: pbar(7),rhobar(7),vxbar(7),vybar(7),vzbar(7)
REAL*8  :: cbar(7),ebar(7)
REAL*8  :: ft,fts,flat

! This setting is compulsory!
indexj = j
flt=0.0
sk2=0.0
! Calculate sj (1)
   if(pbar(indexj+1)-pbar(indexj-1).lt.0.0)then
      sj=1
   else if(pbar(indexj+1)-pbar(indexj-1).eq.0.0)then
      ft=0.0
      goto 100
   else
      sj=-1
   endif
!
   call rflat(ft,indexj,sj,pbar,ebar,rhobar,vxbar,cbar,0,ngrid)
100       continue
#ifdef PIPPO
!
! j=4+sj is the central zone
!
   indexj=indexj+sj
   if(pbar(indexj+1)-pbar(indexj-1).lt.0.0)then
      sj1=1
   else if(pbar(indexj+1)-pbar(indexj-1).eq.0.0)then
      fts=0.0
      goto 50
   else
      sj1=-1
   endif
   call rflat(fts,indexj,sj1,pbar,ebar,rhobar,vxbar,cbar,1,ngrid)
50        continue
! Calculate the flattening parameter (2)
   flat=max(ft,fts)
#endif
flat=ft
!
return
end
!
!******************************************************************
!******************************************************************
!
SUBROUTINE rflat(ftildej,j,sj,pbar,ebar,rhobar,vxbar,cbar,icont,ngrid)
!
USE vector
USE scalar
!
IMPLICIT NONE
INTEGER :: sj,icont,i,j,ngrid
REAL(KIND=8) :: mmmax
REAL(KIND=8) :: essej,nu1,nu2,nu3,WWj
REAL(KIND=8) :: pbar(7),rhobar(7),vxbar(7),cbar(7),ebar(7)
REAL(KIND=8) :: vgrad,pgrad,egrad,omegatd,sigmatd,omega_loc
REAL(KIND=8) :: sigmaj,omegajvgrad,wj,ftildej,wej,kappaj,kappatd,omegaj
REAL(KIND=8) :: sigmatildej,sigma1,sigma2,omegatildej,omega1,omega2,kappatildej,kappa1,kappa2
!
essej = float(sj)
mmmax = 0.0
sigma1 = 0.5d0
sigma2 = 1.0d0
omega1 = 0.75d0
omega2 = 10.0d0
kappa1 = 2.0d0
kappa2 = 0.01d0
nu1=2.00
nu2=0.00
nu3=0.333300
!
! calculate wj (A.1)
!
pgrad=(pbar(j+1)-pbar(j-1))/min(pbar(j+1),pbar(j-1))
vgrad=vxbar(j-1)-vxbar(j+1)
if(pgrad.gt.0.3300.and.vgrad.ge.0.0)then
   wj=1.0
else
   wj=0.0
endif
!
! Equation A.7
!
sigmatildej = abs(pbar(j+1)-pbar(j-1))/min(pbar(j+2),pbar(j-2))
sigmaj = max(mmmax,(sigmatildej-sigma1)/(sigmatildej+sigma2))
!
! Equation A.8
!
if((pbar(j+2)-pbar(j-2)).eq.0.0)then
    pgrad=0.0
else
pgrad=(pbar(j+1)-pbar(j-1))/(pbar(j+2)-pbar(j-2))
endif
if((ebar(j+2)-ebar(j-2)).eq.0.0)then
    egrad=0.0
else
    egrad=(ebar(j+1)-ebar(j-1))/(ebar(j+2)-ebar(j-2))
endif
omegatildej = max(pgrad,egrad)
omegaj = max(mmmax,omega1*(omega2-omegatildej))
!
! Equation A.9 (check tau = 1/rho?)
!
WWj = sqrt( (max(pbar(j+2),pbar(j-2))+0.50*(gamma-1.0)*&
         (pbar(j+2)+pbar(j-2)))*min(rhobar(j+2),rhobar(j-2)) )
wej=essej*WWj/(rhobar(j-2*sj)+vxbar(j-2*sj))

kappatildej = abs((wej-vxbar(j+2*sj)+essej*cbar(j+2*sj))/wej)
kappaj = max(mmmax,(kappatildej-kappa1)/(kappatildej+kappa2))
!
! flattening parameter
!
ftildej = min(omegaj*wj,sigmaj*wj,kappaj)
!
END SUBROUTINE rflat
!
!*****************************************************************
!*****************************************************************
!
#ifdef PIPPO
SUBROUTINE nudet(ngrid)
!
USE vector
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k,ngrid
REAL*8, DIMENSION (:), ALLOCATABLE :: nut
REAL*8 :: nu1,nu2,nu3,du,s,nutmax,nuconf,mmmax
!
! 1-D parameter
!
allocate(nut(ngrid))
mmmax=0.0
nu1=2.00
nu2=0.00
nu3=0.333300
do i=2,ngrid-2
   du=vx(i+2)-vx(i-1)
   nut(i)=nu1*max(mmmax,-du)
enddo
!
! for nut and nu i=i+1/2 !!!
!
do i=3,ngrid-3
   nutmax=max(nut(i-1),nut(i),nut(i+1))
   if(pres(i+1)-pres(i).lt.0.00)then
      s=1.0
   else if(pres(i+1)-pres(i).eq.0.00)then
      s=0.00
   else
      s=-1.00
   endif
   nuconf=(nu2+nu3*max(sk2(i),sk2(i+1)))*dx*s/dt
   nu(i)=s*min(nutmax,nuconf)
enddo
deallocate(nut)
!
END SUBROUTINE nudet
#endif

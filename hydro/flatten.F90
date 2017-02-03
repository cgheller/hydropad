SUBROUTINE flattening(ngrid,nbound)
!
USE vector
USE scalar
!
IMPLICIT NONE
INTEGER :: i,j,k
INTEGER :: ngrid,nbound
INTEGER :: sj,sj1
REAL*8  :: pbar(7),rhobar(7),vxbar(7),vybar(7),vzbar(7)
REAL*8  :: cbar(7),ebar(7)
REAL*8  :: ft,fts

! This setting is compulsory!
flt=0.0
sk2=0.0
! 
! i=1-->boundary
!
i=1+nbound
do j=1,7
      pbar(j)=pres(j)
      rhobar(j)=rho(j)
      vxbar(j)=vx(j)
      vybar(j)=vy(j)
      vzbar(j)=vz(j)
      cbar(j)=sqrt(gamma*pres(j)/rho(j))
      ebar(j)=etot(j)
enddo 
!
! j=4 is the central zone
!
   j=4
   if(pbar(j+1)-pbar(j-1).lt.0.0)then
      sj=1
   else if(pbar(j+1)-pbar(j-1).eq.0.0)then
      flt(i)=1.0
      sk2(i)=0.0
      goto 10
   else
      sj=-1
   endif
!
   call rflat(ft,j,sj,pbar,ebar,rhobar,vxbar,cbar,0,ngrid)
!
! j=4+sj is the central zone
!
   j=j+sj
   if(pbar(j+1)-pbar(j-1).lt.0.0)then
      sj1=1
   else if(pbar(j+1)-pbar(j-1).eq.0.0)then
      fts=1.0
      goto 5
   else
      sj1=-1
   endif
   call rflat(fts,j,sj1,pbar,ebar,rhobar,vxbar,cbar,1,ngrid)
5         continue
   flt(i)=max(ft,fts)
10        continue
!
do i=2+nbound,ngrid-nbound
!
! j=4 is the central zone
!
   do j=1,6
      pbar(j)=pbar(j+1)
      rhobar(j)=rhobar(j+1)
      vxbar(j)=vxbar(j+1)
      vybar(j)=vybar(j+1)
      vzbar(j)=vzbar(j+1)
      cbar(j)=cbar(j+1)
      ebar(j)=ebar(j+1)
   enddo 
      pbar(7)=pres(mod(i+2,ngrid)+1)
      rhobar(7)=rho(mod(i+2,ngrid)+1)
      cbar(7)=sqrt(gamma*pbar(7)/rhobar(7))
      vxbar(7)=vx(mod(i+2,ngrid)+1)
      vybar(7)=vy(mod(i+2,ngrid)+1)
      vzbar(7)=vz(mod(i+2,ngrid)+1)
      ebar(7)=etot(mod(i+2,ngrid)+1)

   j=4
   if(pbar(j+1)-pbar(j-1).lt.0.0)then
      sj=1
   else if(pbar(j+1)-pbar(j-1).eq.0.0)then
      flt(i)=1.0
      sk2(i)=0.0
      goto 100
   else
      sj=-1
   endif
!
   call rflat(ft,j,sj,pbar,ebar,rhobar,vxbar,cbar,0,ngrid)
!
! j=4+sj is the central zone
!
   j=j+sj
   if(pbar(j+1)-pbar(j-1).lt.0.0)then
      sj1=1
   else if(pbar(j+1)-pbar(j-1).eq.0.0)then
      fts=1.0
      goto 50
   else
      sj1=-1
   endif
   call rflat(fts,j,sj1,pbar,ebar,rhobar,vxbar,cbar,1,ngrid)
50        continue
   flt(i)=max(ft,fts)
100       continue
enddo
!
return
end
!
!******************************************************************
!******************************************************************
!
SUBROUTINE rflat(f,j,sj,pbar,ebar,rhobar,vxbar,cbar,icont,ngrid)
!
USE vector
USE scalar
!
IMPLICIT NONE
INTEGER sj,icont,i,j,ngrid
REAL*8 mmmax
REAL*8 pbar(7),rhobar(7),vxbar(7),cbar(7),ebar(7)
REAL*8 pgrad,egrad,omegatd,sigmatd,omega_loc
REAL*8 sigma,vgrad,w,f,wj,wej,s,kappa,kappatd
!
mmmax=0.0
s=float(sj)
pgrad=abs(pbar(j+1)-pbar(j-1))/min(pbar(j+1),pbar(j-1))
vgrad=vxbar(j-1)-vxbar(j+1)
if(pgrad.gt.0.3300.and.vgrad.ge.0.0)then
   w=1.0
else
   w=0.0
endif

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
omegatd=max(pgrad,egrad)
omega_loc=max(mmmax,0.750*(10.0-omegatd))
!
sigmatd=w*abs(pbar(j+2)-pbar(j-2))/min(pbar(j+2),pbar(j-2))
sigma=max(mmmax,(sigmatd-0.50)/(sigmatd+1.0))
!
wj=sqrt( (max(pbar(j+2),pbar(j-2))+0.50*(gamma-1.0)*&
         (pbar(j+2)+pbar(j-2)))*min(rhobar(j+2),rhobar(j-2)) )
wej=s*wj/(rhobar(j-2*sj)+vxbar(j-2*sj))
!
kappatd=abs((wej-vxbar(j+2*sj)+s*cbar(j+2*sj))/wej)
kappa=max(mmmax,(kappatd-2.0)/(kappatd+0.010))
if(icont.eq.0)then
   sk2(i)=sigma*kappa**2
endif
!
f=min(omega_loc*w,sigma*w,kappa)
!
END SUBROUTINE rflat
!
!*****************************************************************
!*****************************************************************
!
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

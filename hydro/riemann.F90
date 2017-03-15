REAL*8 FUNCTION rarwave(f1,f2,vt,vh)
!
!$acc routine seq
IMPLICIT NONE
REAL*8, INTENT(IN) :: f1,f2,vt,vh
!
rarwave=f2+vt*(f1-f2)/(vt-vh)
!
END FUNCTION rarwave
!
!*****************************************************************
!*****************************************************************
!
REAL*8 FUNCTION rhoshock(p2,p1,rho1)
!
USE scalar
!$acc routine seq
!
IMPLICIT NONE
REAL*8, INTENT(IN) :: p2,p1,rho1
!
rhoshock=rho1*(gamma1*p2/p1+1.0)/(gamma1+p2/p1)
!
END FUNCTION rhoshock
!
!*********************************************************************
!*********************************************************************
!
SUBROUTINE riemann(pstar,vstar,pleft,pright,rholeft,rhoright,vxleft,vxright)
!
! This routine uses Van Leer method to solve the Riemann problem
! 
!USE dimension
USE scalar
!$acc routine seq
!
IMPLICIT NONE
REAL(KIND=8)pleft,pright,rholeft,rhoright,vxleft,vxright
INTEGER :: nnewt,nnewtmax
INTEGER :: i,j,k
REAL*8  :: pstar,vstar,vstarl,vstarr,vprimel,vprimer,pstarold
REAL*8  :: err,terr,kappa,cleft,cright,a,b,pcheck
REAL*8  :: pstup
!
! set error parameters
!
err=1.e10
terr=1e-3
nnewt=0
nnewtmax=100
!
! compute first guess for pressure
!
kappa=sqrt(gamma*(pleft+pright)*(rholeft+rhoright)*0.25)
pcheck=0.50*(pleft+pright)+0.50*kappa*(vxleft+vxright)
if(pcheck.gt.0)pstar=pcheck
!
! start newton iteration
!
do while(err.gt.terr)
   nnewt=nnewt+1
   if(nnewt.gt.nnewtmax)then
      write(*,*)'convergenza non raggiunta'
      stop
   endif
   pstarold=pstar
!
   cleft=sqrt(gamma*pleft*rholeft)
   cright=sqrt(gamma*pright*rhoright)
!
   pstup=(1.-(pstar/pleft)**m)
   if(pstar.ge.pleft)then
      a=cleft*sqrt(gf*pstar/pleft+m)
           else if(pstup.eq.0.)then
              a=cleft*sqrt(gf*pstar/pleft+m)
   else
      a=m*cleft*(1.0-pstar/pleft)/(1-(pstar/pleft)**m)
   endif
           pstup=(1.-(pstar/pright)**m)
   if(pstar.ge.pright)then
      b=cright*sqrt(gf*pstar/pright+m)
   else if(pstup.eq.0.)then
      b=cright*sqrt(gf*pstar/pright+m)
   else
      b=m*cright*(1.0-pstar/pright)/(1-(pstar/pright)**m)
   endif
!
   if(pstar.ge.pleft)then
      vprimel=-(a*a+cleft*cleft)/(2*a*a*a)
   else
      vprimel=-(1.0/cleft)*(pstar/pleft)**(-gf)
   endif
   if(pstar.ge.pright)then
      vprimer=(b*b+cright*cright)/(2*b*b*b)
   else
      vprimer=(1.0/cright)*(pstar/pright)**(-gf)
   endif
!
   vstarl=vxleft-(pstar-pleft)/a
   vstarr=vxright+(pstar-pright)/b

   pstar=pstar-(vstarl-vstarr)/(vprimel-vprimer)
   if(pstar.le.0)then
      pstar=pstarold/1000.0
   endif
   vstar=(pleft-pright+a*vxleft+b*vxright)/(a+b)
!
   err=abs(pstarold-pstar)/pstarold
enddo
!
END SUBROUTINE riemann

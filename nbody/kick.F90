SUBROUTINE kick(np,ppos,pvel,ngx,ngy,ngz,gforce)
!
      USE dimension!, ONLY: ngridxpe, ngridype, ngridzpe, nbound
      USE scalar
      USE mpi_inc
!
! local variables
!
IMPLICIT NONE
!
INTEGER :: np,ngx,ngy,ngz
INTEGER :: i1,i2,i3
REAL(KIND=8) :: xi1,xi2,xi3
INTEGER :: i,j,k
INTEGER :: ip1,im1,jp1,j1,jm1,kp1,km1,i3aux,i1aux,i2aux,k1
REAL(KIND=8), INTENT(IN), DIMENSION(3,np) :: ppos
REAL(KIND=8), INTENT(INOUT), DIMENSION(3,np) :: pvel
REAL(KIND=8), INTENT(IN), DIMENSION(3,ngx,ngy,ngz) :: gforce
REAL(KIND=8) :: x1g,x2g,x3g
REAL(KIND=8) :: d1,d2,d3,d4,d5,d6,d7,d8
REAL(KIND=8) :: dd1,dd2,dd3,de1,de2,de3
REAL(KIND=8) :: ffx,ffy,ffz
INTEGER, DIMENSION(ndims) :: partpe
REAL(KIND=8), DIMENSION(ndims) :: boxlaux
INTEGER :: partpeget
!
! move particles using a 3 steps "KDK" method: 
!                 1)  v^{n+1/2}   = v^{n} + h/2 * F^{n}
!                 2)  x^{n+1}     = x^{n} + h v^{n+1/2}
!                 3)  v^{n+1}     = v^{n+1/2} + h/2 * F^{n+1}
! STEP 1+3: half timestep kicks
!
do j=1,np
!
! calculate the cell in which the particle lies in local box coordinates
!
  if(ppos(1,j) == -1)CYCLE
  partpe=0
  partpeget = 0
#ifdef USEMPI
  partpe(1) = int(ppos(1,j)*npesx)
  partpe(2) = int(ppos(2,j)*npesy)
  partpe(3) = int(ppos(3,j)*npesz)
  CALL MPI_CART_RANK(COMM_CART, partpe, partpeget, ierr)
#endif
!
! find the coordinates of a particle within a box in cells units 
!
  boxlaux(1) = real(partpe(1)*ngridxpe)/real(ngridx)
  boxlaux(2) = real(partpe(2)*ngridype)/real(ngridy)
  boxlaux(3) = real(partpe(3)*ngridzpe)/real(ngridz)
  x1g = (ppos(1,j)-boxlaux(1))*real(ngridx)
  x2g = (ppos(2,j)-boxlaux(2))*real(ngridy)
  x3g = (ppos(3,j)-boxlaux(3))*real(ngridz)
!
! find the leading cell for CIC (left cell)
!
  i1 = int(x1g+0.5)+nbound
  j1 = int(x2g+0.5)+nbound
  k1 = int(x3g+0.5)+nbound
!
  ip1=i1+1
  jp1=j1+1
  kp1=k1+1
!
! set the coordinate of the leading cell
!
  if(mype == partpeget)then
!
   xi1=float(i1-nbound)-0.5
   xi2=float(j1-nbound)-0.5
   xi3=float(k1-nbound)-0.5
!
   dd1=x1g-xi1
   dd2=x2g-xi2
   dd3=x3g-xi3
   if(dd1.gt.1.0.or.dd1.lt.0.0)write(*,*)'WARNING: dd1 overflow'!,ppos(1,j),ppos(2,j),ppos(3,j)
   if(dd2.gt.1.0.or.dd2.lt.0.0)write(*,*)'WARNING: dd2 overflow'!,np
   if(dd3.gt.1.0.or.dd3.lt.0.0)write(*,*)'WARNING: dd3 overflow'!,np
   de1=1.0-dd1
   de2=1.0-dd2
   de3=1.0-dd3
   d1= de1*de2*de3
   d2= dd1*de2*de3
   d3= de1*dd2*de3
   d4= dd1*dd2*de3
   d5= de1*de2*dd3
   d6= dd1*de2*dd3
   d7= de1*dd2*dd3
   d8= dd1*dd2*dd3
!
   ffx=gforce(1,i1,j1,k1)*d1+&
      gforce(1,ip1,j1,k1)*d2+&
      gforce(1,i1,jp1,k1)*d3+&
      gforce(1,ip1,jp1,k1)*d4+&
      gforce(1,i1,j1,kp1)*d5+&
      gforce(1,ip1,j1,kp1)*d6+&
      gforce(1,i1,jp1,kp1)*d7+&
      gforce(1,ip1,jp1,kp1)*d8

   ffy=gforce(2,i1,j1,k1)*d1+&
      gforce(2,ip1,j1,k1)*d2+&
      gforce(2,i1,jp1,k1)*d3+&
      gforce(2,ip1,jp1,k1)*d4+&
      gforce(2,i1,j1,kp1)*d5+&
      gforce(2,ip1,j1,kp1)*d6+&
      gforce(2,i1,jp1,kp1)*d7+&
      gforce(2,ip1,jp1,kp1)*d8

   ffz=gforce(3,i1,j1,k1)*d1+&
      gforce(3,ip1,j1,k1)*d2+&
      gforce(3,i1,jp1,k1)*d3+&
      gforce(3,ip1,jp1,k1)*d4+&
      gforce(3,i1,j1,kp1)*d5+&
      gforce(3,ip1,j1,kp1)*d6+&
      gforce(3,i1,jp1,kp1)*d7+&
      gforce(3,ip1,jp1,kp1)*d8
!
! update velocity and position of the particle at t^n+1/2, and set periodic 
! boundary conditions
!
   pvel(1,j)=pvel(1,j)*(1.0-0.5*dat*dtat)+0.5*ffx*dtat
   pvel(2,j)=pvel(2,j)*(1.0-0.5*dat*dtat)+0.5*ffy*dtat
   pvel(3,j)=pvel(3,j)*(1.0-0.5*dat*dtat)+0.5*ffz*dtat
!
  endif

!
! end of loop on particles
!
enddo
!
END SUBROUTINE kick

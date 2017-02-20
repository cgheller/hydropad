!
SUBROUTINE initialize_hydro
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
USE io_mod
#ifdef TRACERS
USE tracers_mod
#endif
!
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
INTEGER :: startx, endx, starty, endy, startz, endz
INTEGER, ALLOCATABLE, DIMENSION(:) :: coordinates
INTEGER :: igrid
INTEGER :: igridx,igridy,igridz
!
t0h0=0.66666666667
!
! initialize time parameters
!
tin=0.0
told=tin
dt=dtinit
!dt=0.00001
dth=0.5*dt
dtold=dt
t=tin+dt
thalf=told+dt/2.0
at=1.0
atnew=at
dat=0.0
ath=1.0
dath=0.0
datnew=0.0
rat=1.0/at
rdtath=dt/ath
omega=1.0
hubble=1.0
redshiftold=0.0
redshift=0.0
thalf=told+dth
!
! initialize a shock tube problem
!

ALLOCATE(coordinates(ndims))
#ifdef USEMPI
CALL MPI_Cart_coords(COMM_CART,mype,ndims,coordinates,ierr)
#endif
rho3d=0.0
p3d=0.0
vx3d=0.0
vy3d=0.0
vz3d=0.0

! Shock tube in the x direction

#ifdef SHOCKX
!$acc data copyin(coordinates) copyout(rho3d,p3d,vx3d,vy3d,vz3d)
startx = coordinates(1)*ngridxpe+1
endx = (coordinates(1)+1)*ngridxpe
!$acc parallel loop collapse(3) gang worker vector
do k=1,nz
 do j=1,ny
  do i=nbound+1,nx-nbound
    igrid = startx+i
    if(igrid < ngridx/2)then
      rho3d(i,j,k) = 1.0d0
      p3d(i,j,k) = 1.0d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    else
      rho3d(i,j,k) = 0.125d0
      p3d(i,j,k) = 0.1d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    endif
  enddo
 enddo
enddo
!$acc end parallel loop
!$acc end data
#endif

#ifdef SHOCKY
startx = coordinates(2)*ngridype+1
endx = (coordinates(2)+1)*ngridype
do i=1,nx
 do k=1,nz
  do j=nbound+1,ny-nbound
    igrid = startx+j
    if(igrid < ngridy/2)then
      rho3d(i,j,k) = 1.0d0
      p3d(i,j,k) = 1.0d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    else
      rho3d(i,j,k) = 0.125d0
      p3d(i,j,k) = 0.1d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    endif
  enddo
 enddo
enddo
#endif


#ifdef SHOCKZ
startx = coordinates(3)*ngridzpe+1
endx = (coordinates(3)+1)*ngridzpe
do i=1,nx
 do j=1,ny
  do k=nbound+1,nz-nbound
    igrid = startx+k
    if(igrid < ngridz/2)then
      rho3d(i,j,k) = 1.0d0
      p3d(i,j,k) = 1.0d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    else
      rho3d(i,j,k) = 0.125d0
      p3d(i,j,k) = 0.1d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    endif
  enddo
 enddo
enddo
#endif

#ifdef SHOCKXYZ
rho3d=real(mype+10)
startx = coordinates(1)*ngridxpe+1
endx = (coordinates(1)+1)*ngridxpe
starty = coordinates(2)*ngridype+1
endy = (coordinates(2)+1)*ngridype
startz = coordinates(3)*ngridzpe+1
endz = (coordinates(3)+1)*ngridzpe
do k=nbound+1,nz-nbound
 do j=nbound+1,ny-nbound
  do i=nbound+1,nx-nbound
    igridx = startx+i
    igridy = starty+j
    igridz = startz+k
    !if(igridx < ngridx/2 .AND. igridy < ngridy/2 .AND. igridz < ngridz/2)then
    !if(igridx+igridy+igridz < ngridx+2)then
    if(igridx+igridy < ngridx+1)then
      rho3d(i,j,k) = 1.0d0
      p3d(i,j,k) = 1.0d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
!    else if(igridx+igridy == ngridx+1)then
!      rho3d(i,j,k) = 0.5625
!      p3d(i,j,k) = 0.55
!      vx3d(i,j,k) = 0.0d0
!      vy3d(i,j,k) = 0.0d0
!      vz3d(i,j,k) = 0.0d0
    else
      rho3d(i,j,k) = 0.125d0
      p3d(i,j,k) = 0.1d0
      vx3d(i,j,k) = 0.0d0
      vy3d(i,j,k) = 0.0d0
      vz3d(i,j,k) = 0.0d0
    endif
  enddo
 enddo
enddo
#endif

#ifdef SEDOV
rho3d = 1.0d0
p3d  = 1.0d0
vx3d = 0.0d0
vy3d = 0.0d0
vz3d = 0.0d0
if(mype .EQ. 0)p3d(nbound+1,nbound+1,nbound+1) = 100000.0
#endif

#ifdef BOUNDARIES
rho3d=-1.0d0
do k=1,nz
 do j=1,ny
  do i=1,nx
    if(i .LE. 2*nbound .AND. i .GT. nbound)rho3d(i,j,k)=1.0+float(npes*mype)
    if(i .GT. nx-2*nbound .AND. i .LE. nx-nbound)rho3d(i,j,k)=2.0+float(npes*mype)
    p3d(i,j,k) = 0.0
    vx3d(i,j,k) = 0.0
    vy3d(i,j,k) = 0.0
    vz3d(i,j,k) = 0.0
  enddo
 enddo
enddo
#endif
#ifdef BOUNDARIES2
rho3d=-1.0
do k=1,nz
 do j=1,ny
  do i=1,nx
    if(j .LE. 2*nbound .AND. j .GT. nbound)rho3d(i,j,k)=3.0+float(npes*mype)
    if(j .GT. ny-2*nbound .AND. j .LE. ny-nbound)rho3d(i,j,k)=4.0+float(npes*mype)
    if(i .LE. 2*nbound .AND. i .GT. nbound)rho3d(i,j,k)=1.0+float(npes*mype)
    if(i .GT. nx-2*nbound .AND. i .LE. nx-nbound)rho3d(i,j,k)=2.0+float(npes*mype)
    p3d(i,j,k) = 0.0
    vx3d(i,j,k) = 0.0
    vy3d(i,j,k) = 0.0
    vz3d(i,j,k) = 0.0
  enddo
 enddo
enddo
#endif
#ifdef BOUNDARIES3
rho3d=-1.0
do k=1,nz
 do j=1,ny
  do i=1,nx
    if(k .LE. 2*nbound .AND. k .GT. nbound)rho3d(i,j,k)=5.0+float(npes*mype)
    if(k .GT. nz-2*nbound .AND. k .LE. nz-nbound)rho3d(i,j,k)=6.0+float(npes*mype)
    if(i .LE. 2*nbound .AND. i .GT. nbound)rho3d(i,j,k)=1.0+float(npes*mype)
    if(i .GT. nx-2*nbound .AND. i .LE. nx-nbound)rho3d(i,j,k)=2.0+float(npes*mype)
    if(j .LE. 2*nbound .AND. j .GT. nbound)rho3d(i,j,k)=3.0+float(npes*mype)
    if(j .GT. ny-2*nbound .AND. j .LE. ny-nbound)rho3d(i,j,k)=4.0+float(npes*mype)
    p3d(i,j,k) = 0.0
    vx3d(i,j,k) = 0.0
    vy3d(i,j,k) = 0.0
    vz3d(i,j,k) = 0.0
  enddo
 enddo
enddo
#endif
#ifdef BOUNDARIES4
rho3d = real(mype)
p3d = 0.0
vx3d = 0.0
vy3d = 0.0
vz3d = 0.0
#endif

END SUBROUTINE initialize_hydro

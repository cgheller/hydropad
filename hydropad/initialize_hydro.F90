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
INTEGER :: igrid
INTEGER :: igridx,igridy,igridz
!
rho3d=0.0
p3d=0.0
vx3d=0.0
vy3d=0.0
vz3d=0.0
!
!$acc data copyin(coordinates) 
!
SELECT CASE (chooseinit)
CASE (0)
   write(*,*)"No initial conditions selected, stopping here..."
   STOP
CASE (1)
! #SHOCKX
!
! Shock tube in the x direction
!
startx = coordinates(1)*ngridxpe+1
endx = (coordinates(1)+1)*ngridxpe
!$acc parallel loop collapse(3) private (startx,endx,igrid) gang worker vector
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

CASE (2)
! #ifdef SHOCKY
startx = coordinates(2)*ngridype+1
endx = (coordinates(2)+1)*ngridype
!$acc parallel loop collapse(3) private (startx,endx,igrid) gang worker vector
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
!$acc end parallel loop

CASE (3)
! #ifdef SHOCKZ
startx = coordinates(3)*ngridzpe+1
endx = (coordinates(3)+1)*ngridzpe
!$acc parallel loop collapse(3) private (startx,endx,igrid) gang worker vector
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
!$acc end parallel loop

CASE (4)
! #ifdef SHOCKXYZ
rho3d=real(mype+10)
startx = coordinates(1)*ngridxpe+1
endx = (coordinates(1)+1)*ngridxpe
starty = coordinates(2)*ngridype+1
endy = (coordinates(2)+1)*ngridype
startz = coordinates(3)*ngridzpe+1
endz = (coordinates(3)+1)*ngridzpe
!$acc parallel loop collapse(3) private (igridx,igridy,igridz) gang worker vector
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
!$acc end parallel loop

CASE (5)
! #ifdef SEDOV
!$acc kernels
rho3d = 1.0d0
p3d  = 1.0d0
vx3d = 0.0d0
vy3d = 0.0d0
vz3d = 0.0d0
if(mype .EQ. 0)p3d(nbound+1,nbound+1,nbound+1) = 100000.0
!$acc end kernels

CASE (6)
!$acc kernels
rho3d = 1.0d0
p3d = rho3d**gamma
vx3d = 0.0d0
vy3d = 0.0d0
vz3d = 0.0d0
!$acc end kernels

CASE (100)
! #ifdef BOUNDARIES
!$acc kernels
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
!$acc end kernels
CASE (101)
! #ifdef BOUNDARIES2
!$acc kernels
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
!$acc end kernels
CASE (102)
! #ifdef BOUNDARIES3
!$acc kernels
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
!$acc end kernels
CASE (103)
! #ifdef BOUNDARIES4
!$acc kernels
rho3d = real(mype)
p3d = 0.0
vx3d = 0.0
vy3d = 0.0
vz3d = 0.0
!$acc end kernels
CASE DEFAULT
   write(*,*)"No initial conditions selected, stopping here..."
   STOP
END SELECT
!
!$acc update host(rho3d,p3d,vx3d,vy3d,vz3d)
!$acc end data

END SUBROUTINE initialize_hydro

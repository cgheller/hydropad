SUBROUTINE alloc_arrays
!
USE dimension
USE matrix
USE vector
USE scalar
#ifdef HEATING
USE heat
#endif
USE mpi_inc
!
IMPLICIT NONE
REAL*8  :: job_size,aux1,aux2,aux3,aux4,aux5
INTEGER :: size_of_real
INTEGER :: error
INTEGER :: nhydrovars
integer :: nxr, nyr, nzr

nxr = int(ngridx/npesx)
nyr = int(ngridy/npesy)
nzr = int(ngridz/npesz)

nxr = int(nxr/ntilex)
nyr = int(nyr/ntiley)
nzr = int(nzr/ntilez)

!nx = nxr + nbound
!ny = nyr + nbound
!nz = nzr + nbound
nx = ngridxpe + nbound + nbound
ny = ngridype + nbound + nbound
nz = ngridzpe + nbound + nbound

!
! allocate baryonic and dark matter basic three dimensional vectors
!
INQUIRE(iolength=size_of_real)job_size
aux1=0.0
aux2=0.0
aux3=0.0

nhydrovars = 10
#ifdef STENCIL
nhydrovars = 15
#endif

aux1=nhydrovars*size_of_real*real(nx*ny*nz)*real(ntilex*ntiley*ntilez)/1e6
if(mype.eq.0)then
  write(*,*)
  write(*,*)'Hydro Size (Mbyte/PE)...',aux1
  write(*,*)'Hydro Size (Mbyte)......',aux1*npes
  write(*,*)
endif

ALLOCATE (p3d(nx,ny,nz),STAT=error)
ALLOCATE (rho3d(nx,ny,nz),STAT=error)
ALLOCATE (vx3d(nx,ny,nz),STAT=error)
ALLOCATE (vy3d(nx,ny,nz),STAT=error)
ALLOCATE (vz3d(nx,ny,nz),STAT=error)
#ifndef STENCIL
ALLOCATE (nes3d(nx,ny,nz),STAT=error)
ALLOCATE (cho3d(nx,ny,nz),STAT=error)
ALLOCATE (cho3dnew(nx,ny,nz),STAT=error)
ALLOCATE (p3dnew(nx,ny,nz),STAT=error)
ALLOCATE (rho3dnew(nx,ny,nz),STAT=error)
ALLOCATE (vx3dnew(nx,ny,nz),STAT=error)
ALLOCATE (vy3dnew(nx,ny,nz),STAT=error)
ALLOCATE (vz3dnew(nx,ny,nz),STAT=error)
ALLOCATE (ttt(nx,ny,nz),STAT=error)
ALLOCATE (pold3d(nx,ny,nz),STAT=error)
ALLOCATE (vxold(nx,ny,nz),STAT=error)
ALLOCATE (vyold(nx,ny,nz),STAT=error)
ALLOCATE (vzold(nx,ny,nz),STAT=error)
ALLOCATE (rhoold(nx,ny,nz),STAT=error)
#endif
!if(nuv.eq.1)ALLOCATE (sncell(nx,ny,nz),STAT=error)
#ifdef GRAVITY
aux3=5.0*job_size*real(nx*ny*(nz+3))/1e6
if(mype.eq.0)then
  write(*,*)
  write(*,*)'Grav Size (Mbyte/PE)...',aux3
  write(*,*)'Grav Size (Mbyte)......',aux3*npes
  write(*,*)
endif
ALLOCATE (phi3d(nx,ny,nz),STAT=error)
ALLOCATE (phiold3d(nx,ny,nz),STAT=error)
ALLOCATE (gxold(nx,ny,0:nz+2),STAT=error)
ALLOCATE (gyold(nx,ny,0:nz+2),STAT=error)
ALLOCATE (gzold(nx,ny,0:nz+2),STAT=error)
#endif
#ifdef NBODY
aux2=7.0*job_size*real(nparmax)/1e6
if(mype.eq.0)then
  write(*,*)
  write(*,*)'N-Body Size (Mbyte/PE)...',aux2
  write(*,*)'N-Body Size (Mbyte)......',aux2*npes
  write(*,*)
endif
ALLOCATE (rhodm3d(nx,ny,nz),STAT=error)
ALLOCATE (x1(nparmax),STAT=error)
ALLOCATE (x2(nparmax),STAT=error)
ALLOCATE (x3(nparmax),STAT=error)
ALLOCATE (v1(nparmax),STAT=error)
ALLOCATE (v2(nparmax),STAT=error)
ALLOCATE (v3(nparmax),11STAT=error)

#endif
#ifdef FFT
ALLOCATE (grav_shap(0:n11,0:n11,0:n11pe),STAT=error)
ALLOCATE (ak(0:n11),STAT=error)
ALLOCATE (akq(0:n11),STAT=error)
#endif
!
if(error.ne.0)then
     write(*,*)mype,' failed in allocating array in ALLOC_ARRAYS'
     stop
endif
!
!$acc kernels
p3d=0.0
rho3d=0.0
vx3d=0.0
vy3d=0.0
vz3d=0.0
!$acc end kernels
#ifndef STENCIL
nes3d=0.0
cho3d=0.0
cho3dnew=0.0
p3dnew=0.0
rho3dnew=0.0
vx3dnew=0.0
vy3dnew=0.0
vz3dnew=0.0
ttt=0.0
pold3d=0.0
vxold=0.0
vyold=0.0
vzold=0.0
rhoold=0.0
#endif
!if(nuv.eq.1)sncell=0.0
#ifdef GRAVITY
phi3d=0.0
phiold3d=0.0
gxold=0.0
gyold=0.0
gzold=0.0
#endif
#ifdef NBODY
rhodm3d=0.0
x1=0.0
x2=0.0
x3=0.0
v1=0.0
v2=0.0
v3=0.0
#endif
#ifdef FFT
grav_shap=0.0
ak=0.0
akq=0.0
#endif
!
END SUBROUTINE alloc_arrays

SUBROUTINE init_sys
!
USE dimension
USE scalar
USE mpi_inc 
!
IMPLICIT NONE
!
! local variables
!
INTEGER :: i
NAMELIST/MESH/ngridx,ngridy,ngridz,nbound
NAMELIST/SETUP/chooseinit
NAMELIST/MPI/npesx,npesy,npesz
NAMELIST/TIME/maxsteps,dtinit,tstop,startdump
NAMELIST/PPM/flatvalue,mingradflat
NAMELIST/NB/npart
NAMELIST/IO/output_pe
NAMELIST/NUMERICS/cour,ca,eta1,eta2,dmax,norm
NAMELIST/PHYSICS/boltz,mh,gamma
NAMELIST/COSMOLOGY/omega_dm,omega_bm,omega_lambda,hnow,boxMpc_over_h,tempnow
!
! Read simulation parameters
!
#ifdef USEMPI
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
do i = 0,npes-1
 if (i .EQ. mype)then
   OPEN(UNIT=500,FILE="input.nml",STATUS="old")
   READ(500,NML=MESH)
   READ(500,NML=SETUP)
#ifdef USEMPI
   READ(500,NML=MPI)
#endif
   READ(500,NML=TIME)
   READ(500,NML=PPM)
   READ(500,NML=NB)
   READ(500,NML=IO)
   READ(500,NML=NUMERICS)
   READ(500,NML=PHYSICS)
   READ(500,NML=COSMOLOGY)

   CLOSE(500)
 endif
#ifdef USEMPI
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
enddo
!
! set-up the domain decomposition
!
ngridxpe = ngridx/npesx
ngridype = ngridy/npesy
ngridzpe = ngridz/npesz
ntilex = 1
ntiley = 1
ntilez = 1
if(npart == -1) npart = ngridx*ngridy*ngridz
npartpe = npart/(npesx*npesy*npesz)
!
! allocate persistent arrays
!
CALL alloc_arrays
!
!
! Create Cartesian MPI framework
!
#ifdef USEMPI
if(npes .ne. npesx*npesy*npesz)then
   if(mype .EQ. 0)write(*,*)"The number of processors does not fit the 3D distribution",npes, npesx, npesy, npesz
   CALL MPI_FINALIZE(ierr)
   stop
endif
#endif
ALLOCATE(dims(ndims))
ALLOCATE(periods(ndims))
ALLOCATE(neighbour(2*ndims))
ALLOCATE(coordinates(ndims))
ALLOCATE(boxl(ndims))
ALLOCATE(boxsize(ndims))

dims(1) = npesx
dims(2) = npesy
dims(3) = npesz
periods(1) = .true.
periods(2) = .true.
periods(3) = .true.
reorder = 0
neighbour = 0
coordinates = 0
boxl = 0.0
boxsize(1) = real(ngridx)/real(max(ngridx,ngridy,ngridz))
boxsize(2) = real(ngridy)/real(max(ngridx,ngridy,ngridz))
boxsize(3) = real(ngridz)/real(max(ngridx,ngridy,ngridz))

#ifdef USEMPI
CALL MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,COMM_CART,ierr)
CALL MPI_Cart_coords(COMM_CART,mype,ndims,coordinates,ierr)

boxl(1) = real(coordinates(1)*ngridxpe)/real(ngridx)
boxl(2) = real(coordinates(2)*ngridype)/real(ngridy)
boxl(3) = real(coordinates(3)*ngridzpe)/real(ngridz)
boxsize(1) = real(ngridxpe)/real(ngridx)
boxsize(2) = real(ngridype)/real(ngridy)
boxsize(3) = real(ngridzpe)/real(ngridz)

if(mype .EQ. output_pe) write(*,*)"MPI Rank ", mype, " has the following neighbours (x:l/r, y:d/u, z:f/r):"
do i=0,ndims-1
    CALL MPI_Cart_shift(COMM_CART,i,1,neighbour(2*i+1),neighbour(2*(i+1)),ierr)
    if(mype .EQ. output_pe) write(*,*)neighbour(2*i+1),neighbour(2*(i+1))
enddo
!
! Create new datatypes (for communication)
! creating x-z planes:
!
CALL MPI_Type_vector(nz, nx*nbound, nx*ny, MPI_DOUBLE_PRECISION, planexznbound, ierr)
CALL MPI_TYPE_COMMIT(planexznbound, ierr)  
!
! creating y-z planes
!
CALL MPI_Type_vector(ny*nz, nbound, nx, MPI_DOUBLE_PRECISION, planeyznbound, ierr)
CALL MPI_TYPE_COMMIT(planeyznbound, ierr)  
!
#else
neighbour = 0
#endif
!
! Now the real initialization stuff
!
CALL initialize_parameters
#ifdef HYDRO
CALL initialize_hydro
#endif
#ifdef NBODY
CALL initialize_nbody
#endif
#ifdef GRAVITY
CALL initialize_gravity
#endif
!
! create ghost regions
!
#ifdef HYDRO
CALL exchange_mesh
#endif
!
nrot=1
!
END SUBROUTINE init_sys

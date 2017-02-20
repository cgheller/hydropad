
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
NAMELIST/MPI/npesx,npesy,npesz
NAMELIST/TIME/maxsteps,dtinit
NAMELIST/HYDRO/flatvalue,mingradflat
NAMELIST/IO/output_pe
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
#ifdef USEMPI
   READ(500,NML=MPI)
#endif
   READ(500,NML=TIME)
   READ(500,NML=HYDRO)
   READ(500,NML=IO)

   CLOSE(500)
 endif
#ifdef USEMPI
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
enddo
ngridxpe = ngridx/npesx
ngridype = ngridy/npesy
ngridzpe = ngridz/npesz
ntilex = 1
ntiley = 1
ntilez = 1
!
! set the geometry, allocate basic arrays
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
ALLOCATE(coords(ndims))
dims(1) = npesx
dims(2) = npesy
dims(3) = npesz
periods(1) = .true.
periods(2) = .true.
periods(3) = .true.
reorder = 0
neighbour = 0

#ifdef USEMPI
CALL MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,COMM_CART,ierr)

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


#endif
!!!CALL indata(startdump)
! CLA: ALL FOLLOWING ROUTINES MUST BE REPLACED
CALL indata_simple
CALL initialize_hydro
CALL exchange_mesh


nrot=1

END SUBROUTINE init_sys

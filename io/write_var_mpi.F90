SUBROUTINE write_var_mpi(auxvar,nx0,ny0,nz0,filename,mype)
!
! update ghost regions
!
USE dimension
USE mpi_inc
IMPLICIT NONE
INTEGER :: nx0,ny0,nz0
REAL(KIND=8), DIMENSION(nx0,ny0,nz0) :: auxvar
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: writevar
INTEGER, ALLOCATABLE, DIMENSION(:) :: sizes
INTEGER, ALLOCATABLE, DIMENSION(:) :: subsizes
INTEGER, ALLOCATABLE, DIMENSION(:) :: starts
!INTEGER, ALLOCATABLE, DIMENSION(:) :: coordinates
INTEGER :: mype,pe
CHARACTER(len=30)filename
!
! local variables
!
INTEGER :: i,j,k,istart
!
INTEGER :: filetype,fhandle,loc_dim
INTEGER(KIND=MPI_OFFSET_KIND)displacement
!
ALLOCATE(sizes(ndims))
ALLOCATE(subsizes(ndims))
ALLOCATE(starts(ndims))
!ALLOCATE(coordinates(ndims))
!
! get task position
!
!CALL MPI_Cart_coords(COMM_CART,mype,ndims,coordinates,ierr)
! 
! Define the file data layout
!
sizes(1) = ngridx
sizes(2) = ngridy
sizes(3) = ngridz
subsizes(1) = ngridxpe
subsizes(2) = ngridype
subsizes(3) = ngridzpe
! starts begin at 0
starts(1) = coordinates(1) * ngridxpe
starts(2) = coordinates(2) * ngridype
starts(3) = coordinates(3) * ngridzpe
loc_dim = ngridxpe*ngridype*ngridzpe
!
! create the datatype
!
CALL MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL, filetype, ierr)
CALL MPI_Type_commit(filetype, ierr)
!
! create the variable to write
!
ALLOCATE(writevar(ngridxpe,ngridype,ngridzpe))
istart=nbound+1
writevar = real(auxvar(istart:istart+ngridxpe-1,istart:istart+ngridype-1,istart:istart+ngridzpe-1))
!
! IO stuff
!
displacement = 0
CALL MPI_File_Open(COMM_CART, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fhandle, ierr)
CALL MPI_File_Set_View(fhandle, displacement, MPI_REAL, filetype, "native", MPI_INFO_NULL, ierr)
CALL MPI_File_Write_All(fhandle, writevar, loc_dim, MPI_REAL, status, ierr)
CALL MPI_File_Close(fhandle, ierr)
!
CALL MPI_Type_free(filetype, ierr)
DEALLOCATE(writevar)
!
END SUBROUTINE write_var_mpi

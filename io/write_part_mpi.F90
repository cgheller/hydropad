SUBROUTINE write_part_mpi(auxvar1,np,filename,mype)
!
! update ghost regions
!
USE dimension
USE mpi_inc
IMPLICIT NONE
INTEGER :: np,ndim0
REAL(KIND=8), DIMENSION(3,np) :: auxvar1
REAL(KIND=4), ALLOCATABLE, DIMENSION(:) :: writevar
INTEGER, ALLOCATABLE, DIMENSION(:) :: sizes
INTEGER, ALLOCATABLE, DIMENSION(:) :: subsizes
INTEGER, ALLOCATABLE, DIMENSION(:) :: starts
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

ndim0=1
ALLOCATE(sizes(ndim0))
ALLOCATE(subsizes(ndim0))
ALLOCATE(starts(ndim0))
! 
! Define the file data layout
!
sizes(1) = 3*npart
subsizes(1) = 3*np
! starts begin at 0
starts(1) = mype*3*np
loc_dim = 3*np
!
! create the datatype
!
CALL MPI_Type_create_subarray(ndim0, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL, filetype, ierr)
CALL MPI_Type_commit(filetype, ierr)
!
! create the variable to write
!
ALLOCATE(writevar(3*np))
istart=1
do i = 1,np
   writevar(j+1) = real(auxvar1(1,i))
   writevar(j+2) = real(auxvar1(2,i))
   writevar(j+3) = real(auxvar1(3,i))
   j=3*(i-1)
enddo
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
END SUBROUTINE write_part_mpi

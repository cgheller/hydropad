#include "hydrompi.def"
!
SUBROUTINE exchange_var(auxvar,nx,ny,nz,nbound,request,nreq,istart,itag)
!
! update ghost regions
!
USE mpi_inc
REAL*8, DIMENSION(nx,ny,nz) :: auxvar
!
! local variables
!
INTEGER :: i,j,k,istart,itag
INTEGER, DIMENSION(nreq) :: commtag
INTEGER, DIMENSION(nreq) :: request
INTEGER :: ytype, ztype
LOGICAL :: smalltiles
!
! Exchange data in the z direction (x-y planes)
!
CALL MPI_Isend(auxvar(1,1,nz-2*nbound+1), nx*ny*nbound, MPI_DOUBLE_PRECISION, neighbour(6), itag, COMM_CART, request(istart), ierr)
istart = istart+1
CALL MPI_Irecv(auxvar(1,1,1),             nx*ny*nbound, MPI_DOUBLE_PRECISION, neighbour(5), itag, COMM_CART, request(istart), ierr)
istart = istart+1
itag = itag+1
CALL MPI_Isend(auxvar(1,1,nbound+1),      nx*ny*nbound, MPI_DOUBLE_PRECISION, neighbour(5), itag, COMM_CART, request(istart), ierr)
istart = istart+1
CALL MPI_Irecv(auxvar(1,1,nz-nbound+1),   nx*ny*nbound, MPI_DOUBLE_PRECISION, neighbour(6), itag, COMM_CART, request(istart), ierr)
istart = istart+1
itag = itag+1
!
! Exchange data in the y direction (x-z planes)
!
CALL MPI_Isend(auxvar(1,ny-2*nbound+1,1), 1, planexznbound, neighbour(4), itag, COMM_CART, request(istart), ierr)
istart = istart+1
CALL MPI_Irecv(auxvar(1,1,1),             1, planexznbound, neighbour(3), itag, COMM_CART, request(istart), ierr)
istart = istart+1
itag = itag+1
CALL MPI_Isend(auxvar(1,nbound+1,1),      1, planexznbound, neighbour(3), itag, COMM_CART, request(istart), ierr)
istart = istart+1
CALL MPI_Irecv(auxvar(1,ny-nbound+1,1),   1, planexznbound, neighbour(4), itag, COMM_CART, request(istart), ierr)
istart = istart+1
itag = itag+1
!
! Exchange data in the x direction (y-z planes)
!
CALL MPI_Isend(auxvar(nx-2*nbound+1,1,1), 1, planeyznbound, neighbour(2), itag, COMM_CART, request(istart), ierr)
istart = istart+1
CALL MPI_Irecv(auxvar(1,1,1),             1, planeyznbound, neighbour(1), itag, COMM_CART, request(istart), ierr)
istart = istart+1
itag = itag+1
CALL MPI_Isend(auxvar(nbound+1,1,1),      1, planeyznbound, neighbour(1), itag, COMM_CART, request(istart), ierr)
istart = istart+1
CALL MPI_Irecv(auxvar(nx-nbound+1,1,1),   1, planeyznbound, neighbour(2), itag, COMM_CART, request(istart), ierr)
istart = istart+1
itag = itag+1

!
END SUBROUTINE exchange_var

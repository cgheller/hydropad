
SUBROUTINE calculate_thermo(tm,pmax0)
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
!
IMPLICIT NONE
!
INTEGER :: i,j,k
INTEGER :: ndata
REAL*8, INTENT(OUT) :: tm,pmax0
REAL*8  :: summa,summa_pe,pmax
!
ndata=nxnynz*npes
summa=0.0
summa_pe=0.0
pmax=0.0
pmax0=0.0
!
do k=1,nz
 do j=1,ny
  do i=1,nx
     ttt(i,j,k)=tfact*p3d(i,j,k)/rho3d(i,j,k)
     summa_pe=summa_pe+ttt(i,j,k)
     pmax=max(pmax,p3d(i,j,k))
  enddo
 enddo
enddo
!
#ifdef USEMPI
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
CALL MPI_Allreduce(summa_pe,summa,1,MPI_DOUBLE_PRECISION,&
                   MPI_SUM,MPI_COMM_WORLD,ierr)
!
CALL MPI_Allreduce(pmax,pmax0,1,MPI_DOUBLE_PRECISION,&
                   MPI_MAX,MPI_COMM_WORLD,ierr)
!
#else
summa=summa_pe
pmax0=pmax
#endif
tm=summa/float(ndata)
!
END SUBROUTINE calculate_thermo	

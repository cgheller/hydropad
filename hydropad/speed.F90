
SUBROUTINE speed
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
!
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
REAL*8 :: vmax,velmaxaux
!
velmax=0.0
velmaxaux=0.0
vmax=0.0
!
#ifdef HYDRO
do k=nbound+1,nz-nbound
   do j=nbound+1,ny-nbound
      do i=nbound+1,nx-nbound
          vsound=sqrt(gamma*p3d(i,j,k)/rho3d(i,j,k))
          vmax=max(vmax,abs(vx3d(i,j,k))+vsound,&
                        abs(vy3d(i,j,k))+vsound,&
                        abs(vz3d(i,j,k))+vsound)
      enddo
   enddo
enddo
!
#ifdef USEMPI
CALL MPI_Allreduce(vmax,velmaxaux,1,MPI_DOUBLE_PRECISION,&
                   MPI_MAX,MPI_COMM_WORLD,ierr)
velmax=velmaxaux
#else
velmax=vmax
#endif
!
#endif
!
vmax = 0.0
#ifdef NBODY
do i = 1,npartpe
  vmax = max(vmax,abs(pvel(1,i)),abs(pvel(2,i)),abs(pvel(3,i)))
enddo
#ifdef USEMPI
CALL MPI_Allreduce(vmax,velmaxaux,1,MPI_DOUBLE_PRECISION,&
                   MPI_MAX,MPI_COMM_WORLD,ierr)
velmax=max(velmax,velmaxaux)
#else
velmax=max(velmax,vmax)
#endif
!
#endif

!
END SUBROUTINE speed

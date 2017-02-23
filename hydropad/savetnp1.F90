SUBROUTINE savetnp1
!
USE dimension
USE matrix
!
IMPLICIT NONE
INTEGER :: i,j,k
!$acc kernels
do k=1,nz
do j=1,nz
do i=1,nz
  rho3dnew(i,j,k) = rho3d(i,j,k)
  vx3dnew(i,j,k) = vx3d(i,j,k)
  vy3dnew(i,j,k) = vy3d(i,j,k)
  vz3dnew(i,j,k) = vz3d(i,j,k)
  p3dnew(i,j,k) = p3d(i,j,k)
  cho3dnew(i,j,k) = cho3d(i,j,k)
enddo
enddo
enddo
!$acc end kernels
!
END SUBROUTINE savetnp1

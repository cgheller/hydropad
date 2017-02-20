SUBROUTINE savetnp1
!
USE dimension
USE matrix
!
IMPLICIT NONE
!$acc kernels
rho3dnew = rho3d
vx3dnew = vx3d
vy3dnew = vy3d
vz3dnew = vz3d
p3dnew = p3d
cho3dnew = cho3d
!$acc end kernels
!
END SUBROUTINE savetnp1

SUBROUTINE savetn
!
USE dimension
USE matrix
!
IMPLICIT NONE
!$acc kernels
rho3d = rho3dnew
vx3d = vx3dnew
vy3d = vy3dnew
vz3d = vz3dnew
p3d = p3dnew
cho3d = cho3dnew
!$acc end kernels
!
END SUBROUTINE savetn

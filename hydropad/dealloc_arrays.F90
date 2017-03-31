SUBROUTINE dealloc_arrays
!
USE dimension
USE matrix
IMPLICIT NONE
!
#ifdef HYDRO
DEALLOCATE (p3d)
DEALLOCATE (rho3d)
DEALLOCATE (vx3d)
DEALLOCATE (vy3d)
DEALLOCATE (vz3d)
#ifndef STENCIL
DEALLOCATE (nes3d)
DEALLOCATE (cho3d)
DEALLOCATE (cho3dnew)
DEALLOCATE (p3dnew)
DEALLOCATE (rho3dnew)
DEALLOCATE (vx3dnew)
DEALLOCATE (vy3dnew)
DEALLOCATE (vz3dnew)
DEALLOCATE (ttt)
DEALLOCATE (pold3d)
DEALLOCATE (vxold)
DEALLOCATE (vyold)
DEALLOCATE (vzold)
DEALLOCATE (rhoold)
#endif
#endif
#ifdef GRAVITY
DEALLOCATE (phi3d)
DEALLOCATE (phiold3d)
DEALLOCATE (gforce)
#endif
#ifdef NBODY
DEALLOCATE (rhodm3d)
DEALLOCATE (ppos)
DEALLOCATE (pvel)
#endif
#ifdef FFT
DEALLOCATE (grav_shap)
#endif
!
END SUBROUTINE dealloc_arrays

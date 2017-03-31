SUBROUTINE dealloc_particles
!
USE nbody_mod
IMPLICIT NONE
!
DEALLOCATE(xrecv)
DEALLOCATE(xsend)
DEALLOCATE(vrecv)
DEALLOCATE(vsend)
DEALLOCATE(nsendpe)
DEALLOCATE(nrecvpe)
DEALLOCATE(indexsend)
DEALLOCATE(psend)
!
END SUBROUTINE dealloc_particles

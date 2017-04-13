SUBROUTINE dealloc_particles
!
USE nbody_mod
IMPLICIT NONE
!
DEALLOCATE(listofparticles)
DEALLOCATE(xrecv)
DEALLOCATE(xsend)
DEALLOCATE(vrecv)
DEALLOCATE(vsend)
DEALLOCATE(nsendpe)
DEALLOCATE(nrecvpe)
DEALLOCATE(indexsend)
DEALLOCATE(psend)
DEALLOCATE(precv)
!
END SUBROUTINE dealloc_particles

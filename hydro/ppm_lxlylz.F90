SUBROUTINE ppm_lxlylz
!
USE matrix
USE scalar
USE mpi_inc
!
IMPLICIT NONE
!
SELECT CASE (nrot)
CASE (1)
  if(mype == 0)write(*,*)"Integration order x-y-z" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
  call savetn
  call exchange_mesh

CASE (2)
  if(mype == 0)write(*,*)"Integration order z-y-x" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
  call savetn
  call exchange_mesh
CASE (3)
  if(mype == 0)write(*,*)"Integration order y-z-x" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
  call savetn
  call exchange_mesh
CASE (4)
  if(mype == 0)write(*,*)"Integration order x-z-y" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
  call savetn
  call exchange_mesh
CASE (5)
  if(mype == 0)write(*,*)"Integration order z-x-y" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
  call savetn
  call exchange_mesh
CASE (6)
  if(mype == 0)write(*,*)"Integration order y-x-z" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
  call savetn
  call exchange_mesh
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
  call savetn
  call exchange_mesh
END SELECT
!
END SUBROUTINE ppm_lxlylz

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
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)

CASE (2)
  if(mype == 0)write(*,*)"Integration order z-y-x" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
CASE (3)
  if(mype == 0)write(*,*)"Integration order y-z-x" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
CASE (4)
  if(mype == 0)write(*,*)"Integration order x-z-y" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
CASE (5)
  if(mype == 0)write(*,*)"Integration order z-x-y" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
CASE (6)
  if(mype == 0)write(*,*)"Integration order y-x-z" 
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 2)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 1)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                 p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                 nx, ny, nz, nbound, 3)
#ifndef STENCIL
  call savetn
#endif
!$acc update host(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
  call exchange_mesh
!!$acc update device(rho3d,vx3d,vy3d,vz3d,p3d,cho3d)
END SELECT
!
END SUBROUTINE ppm_lxlylz

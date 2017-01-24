SUBROUTINE ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                     p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                     nx, ny, nz, nbound, ndirection)
!
! 1D integration in the 
!
USE vector
USE scalar
!
IMPLICIT NONE
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew
!
INTEGER :: nx, ny, nz, ndirection
INTEGER :: nbound
!
! local variables
!
INTEGER :: i,j,k
REAL*8  :: internal,total

SELECT CASE (ndirection)
! Integration in the x direction
CASE (1)
        write(*,*)"Integrating along x"
        ngrid = nx
        call alloc_vectors(ngrid)
        do k=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
            !write(*,*)"------------"
            do i=1,nx
              pres(i)=p3d(i,j,k)
              rho(i)=rho3d(i,j,k)
              vx(i)=vx3d(i,j,k)
              vy(i)=vy3d(i,j,k)
              vz(i)=vz3d(i,j,k)
              eint(i)=internal(i)
              etot(i)=total(i)
              cho(i)=cho3d(i,j,k)
              nes(i)=nes3d(i,j,k)
              c(i)=sqrt(gamma*p3d(i,j,k)/rho3d(i,j,k))
              !write(*,*)rho3d(i,j,k),vx3d(i,j,k)
            enddo
            call flattening(ngrid,nbound)
            call ppm(nbound)
            !write(*,*)"------------"
            do i=1,nx
              p3dnew(i,j,k)=pres(i)
              rho3dnew(i,j,k)=rho(i)
              vx3dnew(i,j,k)=vx(i)
              vy3dnew(i,j,k)=vy(i)
              vz3dnew(i,j,k)=vz(i)
              cho3dnew(i,j,k)=cho(i)
              !write(*,*)rho3dnew(i,j,k),vx3dnew(i,j,k)
            enddo
          enddo
        enddo
        call dealloc_vectors
! Integration in the y direction
CASE (2)
        write(*,*)"Integrating along y"
        ngrid = ny
        call alloc_vectors(ngrid)
        do k=nbound+1,nz-nbound
          do j=nbound+1,nx-nbound
            !write(*,*)"------------"
            do i=1,ny
              pres(i)=p3d(j,i,k)
              rho(i)=rho3d(j,i,k)
              vx(i)=vy3d(j,i,k)
              vy(i)=vz3d(j,i,k)
              vz(i)=vx3d(j,i,k)
              eint(i)=internal(i)
              etot(i)=total(i)
              cho(i)=cho3d(j,i,k)
              nes(i)=nes3d(j,i,k)
              c(i)=sqrt(gamma*p3d(j,i,k)/rho3d(j,i,k))
              !write(*,*)rho3d(i,j,k),vx3d(i,j,k)
            enddo
            call flattening(ngrid,nbound)
            call ppm(nbound)
            !write(*,*)"------------"
            do i=1,ny
              p3dnew(j,i,k)=pres(i)
              rho3dnew(j,i,k)=rho(i)
              vx3dnew(j,i,k)=vz(i)
              vy3dnew(j,i,k)=vx(i)
              vz3dnew(j,i,k)=vy(i)
              cho3dnew(j,i,k)=cho(i)
              !write(*,*)rho3dnew(i,j,k),vx3dnew(i,j,k)
            enddo
          enddo
        enddo
        call dealloc_vectors
! Integration in the z direction
CASE (3)
        write(*,*)"Integrating along z"
        ngrid = nz
        call alloc_vectors(ngrid)
        do k=nbound+1,nx-nbound
          do j=nbound+1,ny-nbound
            do i=1,nz
              pres(i)=p3d(k,j,i)
              rho(i)=rho3d(k,j,i)
              vx(i)=vz3d(k,j,i)
              vy(i)=vy3d(k,j,i)
              vz(i)=vx3d(k,j,i)
              eint(i)=internal(i)
              etot(i)=total(i)
              cho(i)=cho3d(k,j,i)
              nes(i)=nes3d(k,j,i)
              c(i)=sqrt(gamma*p3d(k,j,i)/rho3d(k,j,i))
            enddo
            call flattening(ngrid,nbound)
            call ppm(nbound)
            do i=1,nz
              p3dnew(k,j,i)=pres(i)
              rho3dnew(k,j,i)=rho(i)
              vx3dnew(k,j,i)=vz(i)
              vy3dnew(k,j,i)=vy(i)
              vz3dnew(k,j,i)=vx(i)
              cho3dnew(k,j,i)=cho(i)
            enddo
          enddo
        enddo
        call dealloc_vectors
END SELECT
!
END SUBROUTINE ppmsolver

SUBROUTINE ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                     p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                     nx, ny, nz, nbound, ndirection)
!
! 1D integration  
!

USE vector
USE scalar
USE ppm_mod
!
IMPLICIT NONE
INTEGER :: nx, ny, nz, ndirection,nmin,nmax
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew
INTEGER :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
!
INTEGER :: nbound,ngrid
!
! local variables
!
INTEGER :: i,j,k,ippm,iaux
REAL*8  :: internal,total


SELECT CASE (ndirection)
! Integration in the x direction
CASE (1)
        if(mype == 0)write(*,*)"Integrating along x"
        ngrid = nx
        nmin = nbound
        nmax = ngrid-nbound+1
#ifdef STENCIL
        if(mype == 0)write(*,*)"Using STENCIL"
        ngrid = 7
        nmin = 3
        nmax = 5
#endif
!$OMP PARALLEL
        call alloc_vectors(ngrid)
!$OMP END PARALLEL
!$OMP parallel do collapse(2) default(firstprivate) shared(p3d,rho3d,vx3d,vy3d,vz3d,&
!$OMP          p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew,cho3d,nes3d,cho3dnew) private(i,j,k)
        do k=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
#ifdef STENCIL
            do i=nbound+1,nx-nbound
              do ippm=1,7
                iaux=i+ippm-4 
                pres(ippm)=p3d(iaux,j,k)
                rho(ippm)=rho3d(iaux,j,k)
                vx(ippm)=vx3d(iaux,j,k)
                vy(ippm)=vy3d(iaux,j,k)
                vz(ippm)=vz3d(iaux,j,k)
                eint(ippm)=internal(ippm)
                etot(ippm)=total(ippm)
                cho(ippm)=cho3d(iaux,j,k)
                nes(ippm)=nes3d(iaux,j,k)
                c(ippm)=sqrt(gamma*p3d(iaux,j,k)/rho3d(iaux,j,k))
              enddo
              call ppm(ngrid,nbound,nmin,nmax,i)
              p3dnew(i,j,k)=pres(4)
              rho3dnew(i,j,k)=rho(4)
              vx3dnew(i,j,k)=vx(4)
              vy3dnew(i,j,k)=vy(4)
              vz3dnew(i,j,k)=vz(4)
              cho3dnew(i,j,k)=cho(4)
            enddo
#else
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
            enddo
            call ppm(ngrid,nbound,nmin,nmax,i)
            do i=1,nx
              p3dnew(i,j,k)=pres(i)
              rho3dnew(i,j,k)=rho(i)
              vx3dnew(i,j,k)=vx(i)
              vy3dnew(i,j,k)=vy(i)
              vz3dnew(i,j,k)=vz(i)
              cho3dnew(i,j,k)=cho(i)
            enddo

#endif
          enddo
        enddo
!$OMP end parallel do
        call dealloc_vectors
! Integration in the y direction
CASE (2)
        if(mype == 0)write(*,*)"Integrating along y"
        ngrid = ny
        nmin = nbound
        nmax = ngrid-nbound+1
#ifdef STENCIL
        ngrid = 7
        nmin = 3
        nmax = 5
#endif
!$OMP PARALLEL
        call alloc_vectors(ngrid)
!$OMP END PARALLEL
!$OMP parallel do collapse(2) default(firstprivate) shared(p3d,rho3d,vx3d,vy3d,vz3d,&
!$OMP          p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew,cho3d,nes3d,cho3dnew) private(i,j,k)
        do k=nbound+1,nz-nbound
          do j=nbound+1,nx-nbound
#ifdef STENCIL
            do i=nbound+1,ny-nbound
              do ippm=1,7
                iaux=i+ippm-4
                pres(ippm)=p3d(j,iaux,k)
                rho(ippm)=rho3d(j,iaux,k)
                vx(ippm)=vy3d(j,iaux,k)
                vy(ippm)=vz3d(j,iaux,k)
                vz(ippm)=vx3d(j,iaux,k)
                eint(ippm)=internal(ippm)
                etot(ippm)=total(ippm)
                cho(ippm)=cho3d(j,iaux,k)
                nes(ippm)=nes3d(j,iaux,k)
                c(ippm)=sqrt(gamma*p3d(j,iaux,k)/rho3d(j,iaux,k))
              enddo
              call ppm(ngrid,nbound,nmin,nmax,i)
              p3dnew(j,i,k)=pres(4)
              rho3dnew(j,i,k)=rho(4)
              vx3dnew(j,i,k)=vz(4)
              vy3dnew(j,i,k)=vx(4)
              vz3dnew(j,i,k)=vy(4)
              cho3dnew(j,i,k)=cho(4)
            enddo
#else
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
            enddo
            call ppm(ngrid,nbound,nmin,nmax,i)
            do i=1,ny
              p3dnew(j,i,k)=pres(i)
              rho3dnew(j,i,k)=rho(i)
              vx3dnew(j,i,k)=vz(i)
              vy3dnew(j,i,k)=vx(i)
              vz3dnew(j,i,k)=vy(i)
              cho3dnew(j,i,k)=cho(i)
            enddo
#endif
          enddo
        enddo
!$OMP end parallel do
        call dealloc_vectors
! Integration in the z direction
CASE (3)
        if(mype == 0)write(*,*)"Integrating along z"
        ngrid = nz
        nmin = nbound
        nmax = ngrid-nbound+1
#ifdef STENCIL
        ngrid = 7
        nmin = 3
        nmax = 5
#endif

!$OMP PARALLEL
        call alloc_vectors(ngrid)
!$OMP END PARALLEL
!$OMP parallel do collapse(2) default(firstprivate) shared(p3d,rho3d,vx3d,vy3d,vz3d,&
!$OMP          p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew,cho3d,nes3d,cho3dnew) private(i,j,k)
        do k=nbound+1,nx-nbound
          do j=nbound+1,ny-nbound
#ifdef STENCIL
            do i=nbound+1,nz-nbound
              do ippm=1,7
                iaux=i+ippm-4
                pres(ippm)=p3d(k,j,iaux)
                rho(ippm)=rho3d(k,j,iaux)
                vx(ippm)=vz3d(k,j,iaux)
                vy(ippm)=vy3d(k,j,iaux)
                vz(ippm)=vx3d(k,j,iaux)
                eint(ippm)=internal(ippm)
                etot(ippm)=total(ippm)
                cho(ippm)=cho3d(k,j,iaux)
                nes(ippm)=nes3d(k,j,iaux)
                c(ippm)=sqrt(gamma*p3d(k,j,iaux)/rho3d(k,j,iaux))
              enddo
              call ppm(ngrid,nbound,nmin,nmax,i)
              p3dnew(k,j,i)=pres(4)
              rho3dnew(k,j,i)=rho(4)
              vx3dnew(k,j,i)=vz(4)
              vy3dnew(k,j,i)=vy(4)
              vz3dnew(k,j,i)=vx(4)
              cho3dnew(k,j,i)=cho(4)
            enddo
!STOP
#else
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
            call ppm(ngrid,nbound,nmin,nmax,i)
            do i=1,nz
              p3dnew(k,j,i)=pres(i)
              rho3dnew(k,j,i)=rho(i)
              vx3dnew(k,j,i)=vz(i)
              vy3dnew(k,j,i)=vy(i)
              vz3dnew(k,j,i)=vx(i)
              cho3dnew(k,j,i)=cho(i)
            enddo
!STOP
#endif
          enddo
        enddo
!$OMP end parallel do
        call dealloc_vectors
END SELECT
!
END SUBROUTINE ppmsolver

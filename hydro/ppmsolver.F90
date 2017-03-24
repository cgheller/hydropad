SUBROUTINE ppmsolver(p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d, &
                     p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew, &
                     nx, ny, nz, nbound, ndirection)
!
! 1D integration  
!
#define NUMBLOCKS 100
#define NUMTHREADS 128
USE vector
USE scalar
USE mpi_inc
#ifdef PGIACC
USE cudafor
#endif
!$acc routine (ppm) seq
!$acc routine (internal) seq
!$acc routine (total) seq
!$acc routine (order2interp) seq
!$acc routine (riemannsolver) seq
!$acc routine (integrate) seq
!
IMPLICIT NONE
INTEGER :: nx, ny, nz, ndirection,nmin,nmax
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: p3d, rho3d, vx3d, vy3d, vz3d, cho3d, nes3d
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: p3dnew, rho3dnew, vx3dnew, vy3dnew, vz3dnew, cho3dnew
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: pleft,pright
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: rholeft,rhoright
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: vxleft,vxright
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: vyleft,vyright
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: vzleft,vzright
REAL(KIND=8) :: pm,rhom,vxm,vym,vzm,flux(2),etots,eints
REAL(KIND=8) :: tstart,tend
INTEGER :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
!
INTEGER :: nbound,ngrid
INTEGER :: i,j,k,ippm,iaux
REAL*8  :: internal,total
INTEGER :: istat
#ifdef PGIACC
INTEGER(kind=cuda_count_kind) :: val
!
istat = cudaDeviceGetLimit(val, cudaLimitMallocHeapSize)
write(*,*)"cudaLimitMallocHeapSize = ",val
#endif

!
! allocate left and right arrays for Riemann solver
!
#ifdef STENCIL
allocate(pleft(nx,ny,nz))
allocate(pright(nx,ny,nz))
allocate(rholeft(nx,ny,nz))
allocate(rhoright(nx,ny,nz))
allocate(vxleft(nx,ny,nz))
allocate(vxright(nx,ny,nz))
allocate(vyleft(nx,ny,nz))
allocate(vyright(nx,ny,nz))
allocate(vzleft(nx,ny,nz))
allocate(vzright(nx,ny,nz))

!$acc data present(p3d,rho3d,vx3d,vy3d,vz3d,&
!$acc &            p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew)&
!$acc &            create(pleft,pright,rholeft,rhoright,vxleft,vxright,vyleft,vyright,vzleft,vzright)
!$acc update device(mingradflat,flatvalue,gamma,rat,m,rgamma1,rm,rgamma,dmax,rdx,rdtath,dt,dx,dat,at,&
!$acc &            gamma1,gf,eta1,eta2)
#endif
!
tstart = MPI_WTIME()
SELECT CASE (ndirection)
! Integration in the x direction
CASE (1)
        if(mype == 0)write(*,*)"Integrating along x"
        ngrid = nx
        nmin = nbound
        nmax = ngrid-nbound+1
#ifdef STENCIL
        if(mype == 0)write(*,*)"Using STENCIL"
        ngrid = 6
        nmin = 3
        nmax = 4
#endif


!$OMP PARALLEL
        call alloc_vectors(ngrid)
!$OMP END PARALLEL
!$OMP parallel do collapse(2) default(firstprivate) shared(p3d,rho3d,vx3d,vy3d,vz3d,&
!$OMP          p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew) private(i,j,k)

#ifdef STENCIL
!$acc update device (p3d,rho3d,vx3d,vy3d,vz3d)
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,iaux,ippm,pres,rho,vx,vy,vz,g,c) gang vector
        do k=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
            do i=nbound+1,nx-nbound+1
              do ippm=1,ngrid
                iaux=i+ippm-4 
                pres(ippm)=p3d(iaux,j,k)
                rho(ippm)=rho3d(iaux,j,k)
                vx(ippm)=vx3d(iaux,j,k)
                vy(ippm)=vy3d(iaux,j,k)
                vz(ippm)=vz3d(iaux,j,k)
                g(ippm)=0.0
                c(ippm)=sqrt(gamma*pres(ippm)/rho(ippm))
              enddo
              call order2interp(ngrid,nmin,nmax,nbound,pres,rho,vx,vy,vz,g,c,&
                   pleft(i,j,k),pright(i,j,k),rholeft(i,j,k),rhoright(i,j,k),vxleft(i,j,k),vxright(i,j,k),&
                   vyleft(i,j,k),vyright(i,j,k),vzleft(i,j,k),vzright(i,j,k))
            enddo
          enddo
        enddo
!$acc end parallel
!!!!$acc end parallel loop
!!!!!$acc update device (pleft,pright,rholeft,rhoright,vxleft,vxright,vyleft,vyright,vzleft,vzright)
!
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,pm,rhom,vxm,vym,vzm,rho,vx,ippm,iaux) gang vector
        do k=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
            do i=nbound+1,nx-nbound+1
              do ippm=1,ngrid
                iaux=i+ippm-4
                rho(ippm)=rho3d(iaux,j,k)
                vx(ippm)=vx3d(iaux,j,k)
              enddo
              call riemannsolver(ngrid,rho,vx,pleft(i,j,k),pright(i,j,k),rholeft(i,j,k),rhoright(i,j,k),&
                                 vxleft(i,j,k),vxright(i,j,k),&
                                 vyleft(i,j,k),vyright(i,j,k),vzleft(i,j,k),vzright(i,j,k),&
                                 pm,rhom,vxm,vym,vzm)
              pleft(i,j,k)=pm
              rholeft(i,j,k)=rhom
              vxleft(i,j,k)=vxm
              vyleft(i,j,k)=vym
              vzleft(i,j,k)=vzm               
            enddo
          enddo
        enddo
!$acc end parallel
!!!!!$acc end parallel loop
!
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,etots,eints) gang vector
        do k=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
            do i=nbound+1,nx-nbound
              etots = total(p3d(i,j,k),rho3d(i,j,k),vx3d(i,j,k),vy3d(i,j,k),vz3d(i,j,k))
              eints = internal(p3d(i,j,k))
              call integrate(pleft(i,j,k),pleft(i+1,j,k),rholeft(i,j,k),rholeft(i+1,j,k),&
                             vxleft(i,j,k),vxleft(i+1,j,k),vyleft(i,j,k),vyleft(i+1,j,k),&
                             vzleft(i,j,k),vzleft(i+1,j,k),p3d(i,j,k),rho3d(i,j,k),&
                             vx3d(i,j,k),vy3d(i,j,k),vz3d(i,j,k),etots,eints)
            enddo
          enddo
        enddo
!!!!!$acc end parallel loop
!$acc end parallel
!
#else
        do k=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
            do i=1,nx
              pres(i)=p3d(i,j,k)
              rho(i)=rho3d(i,j,k)
              vx(i)=vx3d(i,j,k)
              vy(i)=vy3d(i,j,k)
              vz(i)=vz3d(i,j,k)
              eint(i)=internal(pres(i))
              etot(i)=total(pres(i),rho(i),vx(i),vy(i),vz(i))
              cho(i)=cho3d(i,j,k)
              nes(i)=nes3d(i,j,k)
              c(i)=sqrt(gamma*p3d(i,j,k)/rho3d(i,j,k))
            enddo
            call ppm(ngrid,nbound,nmin,nmax,i,cho,nes,pres,rho,vx,vy,ghalf,vz,g,c,eint,etot,&
                     rhoold1d,vxold1d,phi1,phi0,flt,nu,sk2)
            do i=1,nx
              p3dnew(i,j,k)=pres(i)
              rho3dnew(i,j,k)=rho(i)
              vx3dnew(i,j,k)=vx(i)
              vy3dnew(i,j,k)=vy(i)
              vz3dnew(i,j,k)=vz(i)
              cho3dnew(i,j,k)=cho(i)
            enddo
          enddo
        enddo
#endif
!$OMP end parallel do
        call dealloc_vectors
!
! Integration in the y direction
!
CASE (2)
!GOTO 10
        if(mype == 0)write(*,*)"Integrating along y"
        ngrid = ny
        nmin = nbound
        nmax = ngrid-nbound+1
#ifdef STENCIL
        ngrid = 6
        nmin = 3
        nmax = 4
#endif
!$OMP PARALLEL
        call alloc_vectors(ngrid)
!$OMP END PARALLEL
!$OMP parallel do collapse(2) default(firstprivate) shared(p3d,rho3d,vx3d,vy3d,vz3d,&
!$OMP          p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew) private(i,j,k)

#ifdef STENCIL
!$acc update device (p3d,rho3d,vx3d,vy3d,vz3d)
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,iaux,ippm,pres,rho,vx,vy,vz,g,c) gang vector
!        do k=nbound+1,nz-nbound
!          do j=nbound+1,nx-nbound
!            do i=nbound+1,ny-nbound+1
        do k=nbound+1,nz-nbound
          do i=nbound+1,ny-nbound+1
            do j=nbound+1,nx-nbound
              do ippm=1,ngrid
                iaux=i+ippm-4
                pres(ippm)=p3d(j,iaux,k)
                rho(ippm)=rho3d(j,iaux,k)
                vx(ippm)=vy3d(j,iaux,k)
                vy(ippm)=vx3d(j,iaux,k)
                vz(ippm)=vz3d(j,iaux,k)
                g(ippm)=0.0
                c(ippm)=sqrt(gamma*pres(ippm)/rho(ippm))
              enddo
              call order2interp(ngrid,nmin,nmax,nbound,pres,rho,vx,vy,vz,g,c,&
                   pleft(j,i,k),pright(j,i,k),rholeft(j,i,k),rhoright(j,i,k),vxleft(j,i,k),vxright(j,i,k),&
                   vyleft(j,i,k),vyright(j,i,k),vzleft(j,i,k),vzright(j,i,k))
            enddo
          enddo
        enddo
!$acc end parallel
!!!!!$acc end parallel loop
!!!!!!!$acc update host(pleft,pright,rholeft,rhoright,vxleft,vxright,vyleft,vyright,vzleft,vzright)
!
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,pm,rhom,vxm,vym,vzm,rho,vx,ippm,iaux) gang vector
!        do k=nbound+1,nz-nbound
!          do j=nbound+1,nx-nbound
!            do i=nbound+1,ny-nbound+1
        do k=nbound+1,nz-nbound
          do i=nbound+1,ny-nbound+1
            do j=nbound+1,nx-nbound
              do ippm=1,ngrid
                iaux=i+ippm-4
                rho(ippm)=rho3d(iaux,j,k)
                vx(ippm)=vy3d(iaux,j,k)
              enddo
              call riemannsolver(ngrid,rho,vx,pleft(j,i,k),pright(j,i,k),rholeft(j,i,k),&
                                 rhoright(j,i,k),vxleft(j,i,k),vxright(j,i,k),&
                                 vyleft(j,i,k),vyright(j,i,k),vzleft(j,i,k),vzright(j,i,k),&
                                 pm,rhom,vxm,vym,vzm)
              pleft(j,i,k)=pm
              rholeft(j,i,k)=rhom
              vxleft(j,i,k)=vxm
              vyleft(j,i,k)=vym
              vzleft(j,i,k)=vzm
            enddo
          enddo
        enddo
!$acc end parallel
!!!!!$acc end parallel loop
!
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,etots,eints) gang vector
!        do k=nbound+1,nz-nbound
!          do j=nbound+1,nx-nbound
!            do i=nbound+1,ny-nbound
        do k=nbound+1,nz-nbound
          do i=nbound+1,ny-nbound
            do j=nbound+1,nx-nbound
              etots = total(p3d(j,i,k),rho3d(j,i,k),vx3d(j,i,k),vy3d(j,i,k),vz3d(j,i,k))
              eints = internal(p3d(j,i,k))
              call integrate(pleft(j,i,k),pleft(j,i+1,k),rholeft(j,i,k),rholeft(j,i+1,k),&
                             vxleft(j,i,k),vxleft(j,i+1,k),vyleft(j,i,k),vyleft(j,i+1,k),&
                             vzleft(j,i,k),vzleft(j,i+1,k),p3d(j,i,k),rho3d(j,i,k),&
                             vy3d(j,i,k),vx3d(j,i,k),vz3d(j,i,k),etots,eints)
            enddo
          enddo
        enddo
!!!!!$acc end parallel loop
!$acc end parallel
!
#else
        do k=nbound+1,nz-nbound
          do j=nbound+1,nx-nbound
            do i=1,ny
              pres(i)=p3d(j,i,k)
              rho(i)=rho3d(j,i,k)
              vx(i)=vy3d(j,i,k)
              vy(i)=vx3d(j,i,k)
              vz(i)=vz3d(j,i,k)
              eint(i)=internal(pres(i))
              etot(i)=total(pres(i),rho(i),vx(i),vy(i),vz(i))
              cho(i)=cho3d(j,i,k)
              nes(i)=nes3d(j,i,k)
              c(i)=sqrt(gamma*p3d(j,i,k)/rho3d(j,i,k))
            enddo
            call ppm(ngrid,nbound,nmin,nmax,i,cho,nes,pres,rho,vx,vy,ghalf,vz,g,c,eint,etot,&
                     rhoold1d,vxold1d,phi1,phi0,flt,nu,sk2)
            do i=1,ny
              p3dnew(j,i,k)=pres(i)
              rho3dnew(j,i,k)=rho(i)
              vx3dnew(j,i,k)=vy(i)
              vy3dnew(j,i,k)=vx(i)
              vz3dnew(j,i,k)=vz(i)
              cho3dnew(j,i,k)=cho(i)
            enddo
          enddo
        enddo
#endif
!$OMP end parallel do
        call dealloc_vectors
!
! Integration in the z direction
!
CASE (3)
!GOTO 10
        if(mype == 0)write(*,*)"Integrating along z"
        ngrid = nz
        nmin = nbound
        nmax = ngrid-nbound+1
#ifdef STENCIL
        ngrid = 6
        nmin = 3
        nmax = 4
#endif
!$OMP PARALLEL
        call alloc_vectors(ngrid)
!$OMP END PARALLEL
!$OMP parallel do collapse(2) default(firstprivate) shared(p3d,rho3d,vx3d,vy3d,vz3d,&
!$OMP          p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew) private(i,j,k)

#ifdef STENCIL
!$acc update device (p3d,rho3d,vx3d,vy3d,vz3d)
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,iaux,ippm,pres,rho,vx,vy,vz,g,c) gang vector
!        do k=nbound+1,nx-nbound
!          do j=nbound+1,ny-nbound
!            do i=nbound+1,nz-nbound+1
        do i=nbound+1,nz-nbound+1
          do j=nbound+1,ny-nbound
            do k=nbound+1,nx-nbound
              do ippm=1,ngrid
                iaux=i+ippm-4
                pres(ippm)=p3d(k,j,iaux)
                rho(ippm)=rho3d(k,j,iaux)
                vx(ippm)=vz3d(k,j,iaux)
                vy(ippm)=vy3d(k,j,iaux)
                vz(ippm)=vx3d(k,j,iaux)
                g(ippm)=0.0
                c(ippm)=sqrt(gamma*pres(ippm)/rho(ippm))
              enddo
              call order2interp(ngrid,nmin,nmax,nbound,pres,rho,vx,vy,vz,g,c,&
                   pleft(k,j,i),pright(k,j,i),rholeft(k,j,i),rhoright(k,j,i),vxleft(k,j,i),vxright(k,j,i),&
                   vyleft(k,j,i),vyright(k,j,i),vzleft(k,j,i),vzright(k,j,i))
            enddo
          enddo
        enddo
!$acc end parallel
!!!!!!$acc end parallel loop
!!!!!$acc update host(pleft,pright,rholeft,rhoright,vxleft,vxright,vyleft,vyright,vzleft,vzright)
!
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,pm,rhom,vxm,vym,vzm,rho,vx,ippm,iaux) gang vector
!        do k=nbound+1,nx-nbound
!          do j=nbound+1,ny-nbound
!            do i=nbound+1,nz-nbound+1
        do i=nbound+1,nz-nbound+1
          do j=nbound+1,ny-nbound
            do k=nbound+1,nx-nbound
              do ippm=1,ngrid
                iaux=i+ippm-4
                rho(ippm)=rho3d(iaux,j,k)
                vx(ippm)=vz3d(iaux,j,k)
              enddo
              call riemannsolver(ngrid,rho,vx,pleft(k,j,i),pright(k,j,i),rholeft(k,j,i),&
                                 rhoright(k,j,i),vxleft(k,j,i),vxright(k,j,i),&
                                 vyleft(k,j,i),vyright(k,j,i),vzleft(k,j,i),vzright(k,j,i),&
                                 pm,rhom,vxm,vym,vzm)
              pleft(k,j,i)=pm
              rholeft(k,j,i)=rhom
              vxleft(k,j,i)=vxm
              vyleft(k,j,i)=vym
              vzleft(k,j,i)=vzm
            enddo
          enddo
        enddo
!$acc end parallel
!!!!!$acc end parallel loop
!
!$acc parallel !num_gangs(NUMBLOCKS) vector_length(NUMTHREADS)
!$acc loop collapse(3) private(i,j,k,etots,eints) gang vector
!        do k=nbound+1,nx-nbound
!          do j=nbound+1,ny-nbound
!            do i=nbound+1,nz-nbound
        do i=nbound+1,nz-nbound
          do j=nbound+1,ny-nbound
            do k=nbound+1,nx-nbound
              etots = total(p3d(k,j,i),rho3d(k,j,i),vx3d(k,j,i),vy3d(k,j,i),vz3d(k,j,i))
              eints = internal(p3d(k,j,i))
              call integrate(pleft(k,j,i),pleft(k,j,i+1),rholeft(k,j,i),rholeft(k,j,i+1),&
                             vxleft(k,j,i),vxleft(k,j,i+1),vyleft(k,j,i),vyleft(k,j,i+1),&
                             vzleft(k,j,i),vzleft(k,j,i+1),p3d(k,j,i),rho3d(k,j,i),&
                             vz3d(k,j,i),vy3d(k,j,i),vx3d(k,j,i),etots,eints)
            enddo
          enddo
        enddo
!!!!!$acc end parallel loop
!$acc end parallel
!
#else
        do k=nbound+1,nx-nbound
          do j=nbound+1,ny-nbound
            do i=1,nz
              pres(i)=p3d(k,j,i)
              rho(i)=rho3d(k,j,i)
              vx(i)=vz3d(k,j,i)
              vy(i)=vy3d(k,j,i)
              vz(i)=vx3d(k,j,i)
              eint(i)=internal(pres(i))
              etot(i)=total(pres(i),rho(i),vx(i),vy(i),vz(i))
              cho(i)=cho3d(k,j,i)
              nes(i)=nes3d(k,j,i)
              c(i)=sqrt(gamma*p3d(k,j,i)/rho3d(k,j,i))
            enddo
            call ppm(ngrid,nbound,nmin,nmax,i,cho,nes,pres,rho,vx,vy,ghalf,vz,g,c,eint,etot,&
                     rhoold1d,vxold1d,phi1,phi0,flt,nu,sk2)
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
#endif
!$OMP end parallel do
        call dealloc_vectors
END SELECT
10 CONTINUE
#ifdef STENCIL
!$acc end data
deallocate(pleft)
deallocate(pright)
deallocate(rholeft)
deallocate(rhoright)
deallocate(vxleft)
deallocate(vxright)
deallocate(vyleft)
deallocate(vyright)
deallocate(vzleft)
deallocate(vzright)
#endif
tend = MPI_WTIME()
if(mype == 0)write(*,*)"Elapsed time = ", tend-tstart
!
END SUBROUTINE ppmsolver

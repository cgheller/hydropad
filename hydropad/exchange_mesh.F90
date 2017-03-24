!
SUBROUTINE exchange_mesh
!
! update ghost regions
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
!
! local variables
!
#define NEXCHVAR 100
INTEGER :: i,j,k,istart,itot,itag
INTEGER, DIMENSION(NEXCHVAR) :: request
INTEGER, DIMENSION(MPI_STATUS_SIZE,NEXCHVAR) :: statcomm
!
!
! Exchange data to build ghost regions
!
istart=1
#ifdef USEMPI
itag=10
!variablefilename="rho1.bin"
!if(nstep .eq. 9)call write_var(rho3d,nx,ny,nz,variablefilename,mype,output_pe)
!call write_var_mpi(rho3d,nx,ny,nz,variablefilename,mype)
CALL exchange_var(rho3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(p3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(vx3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(vy3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(vz3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
#ifndef STENCIL
CALL exchange_var(cho3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(nes3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
#endif
itot=istart-1
!
! complete data exchange
!
CALL MPI_Waitall(itot,request,statcomm,ierr)
!
! alternative data exchange if no MPI is used
#else
!
else
do k=1,nbound
 do j=1,ny
  do i=1,nx
    rho3d(i,j,k) = rho3d(i,j,nz-2*nbound+k) 
    p3d(i,j,k) = p3d(i,j,nz-2*nbound+k) 
    vx3d(i,j,k) = vx3d(i,j,nz-2*nbound+k) 
    vy3d(i,j,k) = vy3d(i,j,nz-2*nbound+k) 
    vz3d(i,j,k) = vz3d(i,j,nz-2*nbound+k) 
#ifndef STENCIL
    cho3d(i,j,k) = cho3d(i,j,nz-2*nbound+k) 
    nes3d(i,j,k) = nes3d(i,j,nz-2*nbound+k) 
#endif

    rho3d(i,j,nz-nbound+k) = rho3d(i,j,k+nbound)
    p3d(i,j,nz-nbound+k) = p3d(i,j,k+nbound)
    vx3d(i,j,nz-nbound+k) = vx3d(i,j,k+nbound)
    vy3d(i,j,nz-nbound+k) = vy3d(i,j,k+nbound)
    vz3d(i,j,nz-nbound+k) = vz3d(i,j,k+nbound)
#ifndef STENCIL
    cho3d(i,j,nz-nbound+k) = cho3d(i,j,k+nbound)
    nes3d(i,j,nz-nbound+k) = nes3d(i,j,k+nbound)
#endif
  enddo
 enddo
enddo

do i=1,nbound
 do j=1,ny
  do k=1,nz
    rho3d(i,j,k) = rho3d(nx-2*nbound+i,j,k)
    p3d(i,j,k) = p3d(nx-2*nbound+i,j,k)
    vx3d(i,j,k) = vx3d(nx-2*nbound+i,j,k)
    vy3d(i,j,k) = vy3d(nx-2*nbound+i,j,k)
    vz3d(i,j,k) = vz3d(nx-2*nbound+i,j,k)
#ifndef STENCIL
    cho3d(i,j,k) = cho3d(nx-2*nbound+i,j,k)
    nes3d(i,j,k) = cho3d(nx-2*nbound+i,j,k)
#endif

    rho3d(nx-nbound+i,j,k) = rho3d(i+nbound,j,k)
    p3d(nx-nbound+i,j,k) = p3d(i+nbound,j,k)
    vx3d(nx-nbound+i,j,k) = vx3d(i+nbound,j,k)
    vy3d(nx-nbound+i,j,k) = vy3d(i+nbound,j,k)
    vz3d(nx-nbound+i,j,k) = vz3d(i+nbound,j,k)
#ifndef STENCIL
    cho3d(nx-nbound+i,j,k) = cho3d(i+nbound,j,k)
    nes3d(nx-nbound+i,j,k) = nes3d(i+nbound,j,k)
#endif
  enddo
 enddo
enddo

do j=1,nbound
 do k=1,nz
  do i=1,nx
    rho3d(i,j,k) = rho3d(i,ny-2*nbound+j,k)
    p3d(i,j,k) = p3d(i,ny-2*nbound+j,k)
    vx3d(i,j,k) = vx3d(i,ny-2*nbound+j,k)
    vy3d(i,j,k) = vy3d(i,ny-2*nbound+j,k)
    vz3d(i,j,k) = vz3d(i,ny-2*nbound+j,k)
#ifndef STENCIL
    cho3d(i,j,k) = cho3d(i,ny-2*nbound+j,k)
    nes3d(i,j,k) = nes3d(i,ny-2*nbound+j,k)
#endif

    rho3d(i,ny-nbound+j,k) = rho3d(i,j+nbound,k)
    p3d(i,ny-nbound+j,k) = p3d(i,j+nbound,k)
    vx3d(i,ny-nbound+j,k) = vx3d(i,j+nbound,k)
    vy3d(i,ny-nbound+j,k) = vy3d(i,j+nbound,k)
    vz3d(i,ny-nbound+j,k) = vz3d(i,j+nbound,k)
#ifndef STENCIL
    cho3d(i,ny-nbound+j,k) = cho3d(i,j+nbound,k)
    nes3d(i,ny-nbound+j,k) = nes3d(i,j+nbound,k)
#endif
  enddo
 enddo
enddo
#endif 

END SUBROUTINE exchange_mesh

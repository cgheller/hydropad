SUBROUTINE initialize_gravity
!
USE dimension
USE matrix
USE scalar
USE mpi_inc
!
IMPLICIT NONE
!
! local variables
!
#define NEXCHVAR 100
INTEGER :: i,j,k,istart,itot,itag
INTEGER, DIMENSION(NEXCHVAR) :: request
INTEGER, DIMENSION(MPI_STATUS_SIZE,NEXCHVAR) :: statcomm
REAL(KIND=8) :: deltax,deltay,deltaz,alphagauss,sigma
REAL(KIND=8), DIMENSION(3) :: gcenter
!
alphagauss = 1.0
sigma = 0.3 
gcenter(1) = 0.5
gcenter(2) = 0.5
gcenter(3) = 0.5
!
do k=1,ngridz
 do j=1,ngridy
  do i=1,ngridx
    deltax = (real(i+mype*ngridx)-0.5)*dx-gcenter(1)   
    deltay = (real(j+mype*ngridy)-0.5)*dx-gcenter(2)   
    deltaz = (real(k+mype*ngridz)-0.5)*dx-gcenter(3)   
    phi3d(i+nbound,j+nbound,k+nbound) = alphagauss * exp(-(deltax**2+deltay**2+deltaz**2)/(2.0*sigma**2))
  enddo
 enddo
enddo
!
istart=1
#ifdef USEMPI
itag=10
CALL exchange_var(phi3d,nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
itot=istart-1
CALL MPI_Waitall(itot,request,statcomm,ierr)
#endif
!
do k=nbound+1,ngridz+nbound
 do j=nbound+1,ngridy+nbound
  do i=nbound+1,ngridx+nbound
    gforce(1,i,j,k) = phi3d(i+1,j,k)-phi3d(i-1,j,k)
    gforce(2,i,j,k) = phi3d(i,j+1,k)-phi3d(i,j-1,k)
    gforce(3,i,j,k) = phi3d(i,j,k+1)-phi3d(i,j,k-1)
  enddo
 enddo
enddo
!
istart=1
#ifdef USEMPI
CALL exchange_var(gforce(1,:,:,:),nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(gforce(2,:,:,:),nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
CALL exchange_var(gforce(3,:,:,:),nx,ny,nz,nbound,request,NEXCHVAR,istart,itag)
itot=istart-1
CALL MPI_Waitall(itot,request,statcomm,ierr)
#endif
!
END SUBROUTINE initialize_gravity

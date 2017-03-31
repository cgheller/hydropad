SUBROUTINE initialize_nbody
!
USE dimension
USE matrix
USE scalar
USE mpi_inc
USE io_mod
!
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
INTEGER :: startx, endx, starty, endy, startz, endz
INTEGER :: igrid
INTEGER :: igridx,igridy,igridz
!

do k=1,ngridz
 do j=1,ngridy
  do i=1,ngridx
!
    igrid = i+ngridx*(j-1)+(ngridx*ngridy)*(k-1) + mype*npartpe
!    
    ppos(1,igrid) = ((i-1)+mype*ngridx+0.5)*dx  
    ppos(2,igrid) = ((j-1)+mype*ngridy+0.5)*dx  
    ppos(3,igrid) = ((k-1)+mype*ngridz+0.5)*dx  
    pvel(1,igrid) = 0.0
    pvel(2,igrid) = 0.0
    pvel(3,igrid) = 0.0
    rhodm3d(i,j,k) = amass
  enddo
 enddo
enddo


END SUBROUTINE initialize_nbody

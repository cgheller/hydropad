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

do k=1,ngridzpe
 do j=1,ngridype
  do i=1,ngridxpe
!
    igrid = i+ngridxpe*(j-1)+(ngridxpe*ngridype)*(k-1) 
!    
    ppos(1,igrid) = ((i-1)+coordinates(1)*ngridxpe+0.5)*dx  
    ppos(2,igrid) = ((j-1)+coordinates(2)*ngridype+0.5)*dx  
    ppos(3,igrid) = ((k-1)+coordinates(3)*ngridzpe+0.5)*dx  
    pvel(1,igrid) = 0.0
    pvel(2,igrid) = 0.0
    pvel(3,igrid) = 0.0
    rhodm3d(i,j,k) = amass
  enddo
 enddo
enddo


END SUBROUTINE initialize_nbody

SUBROUTINE write_var(auxvar,nx,ny,nz,filename,mype,pe)
!
! update ghost regions
!
USE mpi_inc
REAL*8, DIMENSION(nx,ny,nz) :: auxvar
INTEGER :: mype,pe
CHARACTER(len=30)filename
!
! local variables
!
INTEGER :: i,j,k
!
if(pe .EQ. mype)then
  open(unit=501,file=filename, access="direct", recl=4*nx*ny*nz)
  write(501,rec=1)real(auxvar)
  close(501)
endif
!
END SUBROUTINE write_var

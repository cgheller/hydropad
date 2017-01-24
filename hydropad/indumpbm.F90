#include "hydrompi.def"

SUBROUTINE indumpbm
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
	USE heat
        USE mpi_inc
!
	IMPLICIT NONE
!
! local variables
!
        CHARACTER :: ch_i1*1, ch_i2*1, ch_i3*1, filename*20
        INTEGER :: nfile,nfilemin,nfilemax,ifile
        INTEGER :: i1,i2,i3,i2aux,iunit
        INTEGER :: i,j,k
!
	k=0
#ifdef DEBUG
	write(*,*)'loading data...',mype
#endif
        nfile=nz
        nfilemin=mype*nfile
        nfilemax=nfilemin+nfile-1
        do ifile=nfilemin,nfilemax
! opening files
          i1=ifile/100
          i2aux=ifile-i1*100
          i2=i2aux/10
          i3=i2aux-i2*10
          i1=i1+48
          i2=i2+48
          i3=i3+48
          ch_i1=char(i1)
          ch_i2=char(i2)
          ch_i3=char(i3)
          filename='outdumpbm.'//ch_i1//ch_i2//ch_i3
          iunit=200+ifile
          open (iunit,file=filename,status='unknown',form='unformatted')
!
	  read(iunit)npl,nrot,nstep
	  read(iunit)at,atnew
	  read(iunit)told,t,dt,dtold
	  read(iunit)ntotsn
	  read(iunit)noutput 
	  k=k+1
          read(iunit)((p3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((rho3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((vx3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((vy3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((vz3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((cho3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((phi3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((phiold3d(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((gxold(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((gyold(i,j,k),i=1,nx),j=1,ny)
          read(iunit)((gzold(i,j,k),i=1,nx),j=1,ny)
	  if(nuv.eq.1)&
          read(iunit)((sncell(i,j,k),i=1,nx),j=1,ny)
!
	  close(iunit)
	enddo
!
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
END SUBROUTINE indumpbm

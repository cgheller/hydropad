!
SUBROUTINE indumpdm
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
!
	IMPLICIT NONE
!
! local variables
!
        character :: ch_i1*1, ch_i2*1, ch_i3*1, filename*20
        integer :: nfile,nfilemin,nfilemax,ifile
        integer :: i1,i2,i3,i2aux,iunit
        integer :: npxpe,iimin,iimax,ics
        INTEGER :: i,j,k
!
        npxpe=nparmax/nz
        nfile=nz
        nfilemin=mype*nfile
        nfilemax=nfilemin+nfile-1
	ics=0
        do ifile=nfilemin,nfilemax
	  iimin=ics*npxpe+1
	  iimax=(ics+1)*npxpe
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
          filename='outdumpdm.'//ch_i1//ch_i2//ch_i3
          iunit=200+ifile
          open(iunit,file=filename,status='unknown',form='unformatted')
!
          read(iunit)(x1(i),i=iimin,iimax)
          read(iunit)(x2(i),i=iimin,iimax)
          read(iunit)(x3(i),i=iimin,iimax)
          read(iunit)(v1(i),i=iimin,iimax)
          read(iunit)(v2(i),i=iimin,iimax)
          read(iunit)(v3(i),i=iimin,iimax)
!
	  ics=ics+1
          close(iunit)
	enddo
!
#ifdef DEBUG
	write(*,*)'loading data...completed',mype
#endif
!
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
END SUBROUTINE indumpdm

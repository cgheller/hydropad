#include "hydrompi.def"
!
SUBROUTINE checkout
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
        USE io_mod
!
! local variables
!
	IMPLICIT NONE
	REAL(kind=4) :: vaux1,vaux2,vaux3
        CHARACTER*1  :: ch1,ch2,ch3,ch4,ch5
        CHARACTER*30 :: filename1,filename6
        INTEGER      :: rec,ii,nlrec1,nrec
        INTEGER      :: iunit1,ijump,njump
        INTEGER      :: i,j,k
	REAL*8       :: xnchk,xnpchk,xnstep,xpart
!
	iunit1=901
	nlrec1=3*nlrec
!
        call getindex10000(nstep,ch1,ch2,ch3,ch4,ch5)
        filename1='check_evol.'//ch1//ch2//ch3//ch4//ch5//'.pff'
!
        do ii=0,npes-1
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
        if(mype.eq.ii)then
!	
        open(unit=iunit1,file=filename1,access='direct',recl=nlrec1)
!
	vaux1=real(nstep)
	vaux2=real(redshift)
        xnpchk=dble(npchk)
	njump=nparmax/npchk
	vaux3=real(npes*njump)
!	write(iunit1,rec=1)vaux1,vaux2,vaux3
	ijump=0
	do i=1,nparmax,npchk
	  ijump=ijump+1
	  nrec=ijump+mype*njump
	  vaux1=real(x1(i))
	  vaux2=real(x2(i))
	  vaux3=real(x3(i))
	  write(iunit1,rec=nrec)vaux1,vaux2,vaux3
	enddo
!
	close(iunit1)
	endif
	enddo
!
END SUBROUTINE checkout

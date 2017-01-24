SUBROUTINE readtables
!
      USE dimension
      USE scalar
      USE ccool
      IMPLICIT NONE

      integer naux,i,j,k,iunit

      iunit=1000+mype
      open(unit=iunit,file='lambdarho.fp',form='unformatted',status='old')

      read(iunit)naux	
      do i=1,rhotable
         read(iunit)trtable(i),mrtable(i),qrtable(i)
      enddo
      close(iunit)

      open(unit=iunit,file='lambdarho2.fp',form='unformatted',status='old')
      read(iunit)naux	
      do i=1,rho2table
         read(iunit)tr2table(i),mr2table(i),qr2table(i)
      enddo
      close(iunit)

      open(unit=iunit,file='lambdarho3.fp',form='unformatted',status='old')
      read(iunit)naux	
      do i=1,rho3table
         read(iunit)tr3table(i),mr3table(i),qr3table(i)
      enddo
      close(iunit)	

END SUBROUTINE readtables 

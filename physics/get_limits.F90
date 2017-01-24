  SUBROUTINE get_limits(flag_telescope)
! returns limits of energy band

  USE scalar
  USE xproj
  IMPLICIT NONE

  real*8 e1,e2,emin,emax,zz
  integer nch,jbin1,jbin2,ich
  
  integer flag_telescope
!
  zz=redshift
  IF (flag_telescope.eq.0) THEN
    e1=0.1
    e2=20
  ELSE IF (flag_telescope.eq.1) THEN
!...CXO ACIS-S                     De=150 eV
    nch=60
    emin=0.5
    emax=9.5
    e1=emin+(emax-emin)/nch*(ich-1)
    e2=emin+(emax-emin)/nch*ich
  ELSE IF (flag_telescope.eq.-1) THEN
!...dummy CXO ACIS-S   
    nch=6
    emin=0.5
    emax=9.5
    e1=emin+(emax-emin)/nch*(ich-1)
    e2=emin+(emax-emin)/nch*ich
  ELSE IF (flag_telescope.eq.2) THEN
!...ASCA SIS                       De/e=0.02*(5.9/e)^1/2
    emin=3.
    emax=12.
    CALL sisband(emin,ich,e1,e2)
  ENDIF
!
  if(mype.eq.0)PRINT*,'Channel energy band lmts :',e1,e2,'keV'
!...spectral bin indices corresponding 
!...to the energy band limits (e1,e2).
  jbin1=INT((e1-powmin)/de)+1
  jbin2=INT((e2-powmin)/de)+1
  if(mype.eq.0)PRINT*,'Energy bins              :',jbin2-jbin1+1
!
!...energy band at the source 
  ez1=e1*(1.+zz)
  ez2=e2*(1.+zz)
  if(mype.eq.0)PRINT*,'Emitted energy band lmts :',ez1,ez2,'keV'
!...spectral bin indices corresponding 
!...to the energy band limits (ez1,ez2) in emission.
  ibin1=INT((ez1-powmin)/de)+1
  ibin2=INT((ez2-powmin)/de)+1
  if(mype.eq.0)PRINT*,'Energy bins              :',ibin2-ibin1+1
!
  RETURN
  END SUBROUTINE get_limits
!..........................................................................
! evaluation of band of the channel ich of the ASCA SIS
 SUBROUTINE sisband(emin,ich,e1,e2)
!      
 IMPLICIT NONE

 INTEGER :: i,ich
 REAL*8  :: c1,c2,aa,e2,emin,e1

 c1=0.02
 c2=5.9
 aa=c1*c1*c2/4.
!      
 e2=emin
 DO i=1,ich
   e1=e2
   e2=e1+aa+sqrt(2*e1*aa+aa*aa)
 ENDDO

 END SUBROUTINE sisband


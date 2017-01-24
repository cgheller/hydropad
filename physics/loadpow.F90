  SUBROUTINE loadpow
!
  USE xproj
  IMPLICIT NONE
!
  INTEGER :: itemp,ichan
  INTEGER :: lblnk
  REAL*8 emean,ebin
!
  CHARACTER(LEN=90) name0
  CHARACTER(LEN=3) suff
!
  DO itemp=0,(ntemp-1)
    CALL itoa(itemp,suff)
    name0=filepow(1:lblnk(filepow))//suff
    OPEN(10,file=name0,status='old')
    DO ichan=1,nchan
      READ(10,*) emean,ebin,pow(itemp,ichan)
    ENDDO
    CLOSE(10)
  ENDDO
  END SUBROUTINE loadpow
!..........................................................................
!
  SUBROUTINE itoa(num,suff)
  CHARACTER(LEN=10) dumfile
  CHARACTER(LEN=3) suff
  integer num
  IF (num.lt.10) THEN
    WRITE(dumfile,'(i1)') num
  ELSE IF (num.lt.100) THEN
    WRITE(dumfile,'(i2)') num
  ELSE 
    WRITE(dumfile,'(i3)') num
  ENDIF
  READ(dumfile,'(a)') suff
  END SUBROUTINE itoa
!..........................................................................
  INTEGER FUNCTION lblnk(char)
!
  CHARACTER (LEN=*) char
  lblnk=index(char,' ')-1
  END FUNCTION lblnk


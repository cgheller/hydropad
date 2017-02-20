REAL*8 FUNCTION internal(pres)
!
! calculate internal energy density
!
	USE dimension
	USE scalar
!$acc routine seq
!
! local variables
!
	IMPLICIT NONE
        REAL(KIND=8)pres
!
	internal=pres*rgamma1
!
END FUNCTION internal

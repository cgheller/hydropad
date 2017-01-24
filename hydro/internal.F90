REAL*8 FUNCTION internal(nnn)
!
! calculate internal energy density
!
	USE dimension
	USE vector
	USE scalar
!
! local variables
!
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nnn
!
	internal=pres(nnn)*rgamma1
!
END FUNCTION internal

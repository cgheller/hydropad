REAL*8 FUNCTION total(nnn)
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
	total=pres(nnn)*rgamma1+0.5*rho(nnn)*(vx(nnn)*vx(nnn)+&
              vy(nnn)*vy(nnn)+vz(nnn)*vz(nnn))
!
END FUNCTION total

SUBROUTINE flux(rneh,T,bcell)
!
	USE xproj
	IMPLICIT NONE
	REAL*8 :: sp_e(ibin1:ibin2)
	REAL*8 :: rneh,T
	REAL*8 :: bcell
!
! select the right spectrum
!
	if(T.ge.Tmin)then
	  CALL get_sp_emiss(rneh,T,sp_e)
!
! integrate over the energy band
!
	  CALL luminosity(sp_e,bcell)
	else
	  bcell=0.0
	endif
!
END SUBROUTINE flux

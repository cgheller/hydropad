  SUBROUTINE luminosity(sp_e,sum)
!...returns luminosity in the energy band of the channel

  USE xproj
  IMPLICIT NONE

  REAL*8 :: sp_e(ibin1:ibin2)
  REAL*8 :: sum,ebin,de_,elow1,elow2
  REAL*8, PARAMETER :: uerg=1.6022e-9
  REAL*8 :: get_de
  INTEGER :: ibin
!
!...integration of differential emissivity on energy and pixels
  sum = 0.
  DO ibin=ibin1,ibin2  
    ebin = powmin + de * (ibin-0.5)           ! channel related energy (keV)
    de_ = get_de(ibin,ibin1,ibin2,ez1,ez2)
    sum= sum + sp_e(ibin) * ebin * uerg * de_ ! emissivity in erg/s/cm^3
  ENDDO
!
  END SUBROUTINE luminosity
!............................................................................
  REAL*8 FUNCTION get_de(ibin,imin,imax,emin,emax)
!...returns energy interval for integration

  USE xproj
  IMPLICIT NONE

  INTEGER ibin,imin,imax
  REAL*8 emin,emax
  real*8 elow1,elow2
!
  elow1=powmin+de*(imin-1)
  elow2=powmin+de*(imax-1)
!
  IF(ibin.eq.imin) THEN
    IF(imin.eq.imax) THEN
      get_de = (emax-emin)
    ELSE
      get_de = (elow1+de-emin)
    ENDIF
  ELSE IF (ibin.lt.imax) THEN
      get_de = de
    ELSE
      get_de = (emax-elow2)
  ENDIF
  END FUNCTION get_de


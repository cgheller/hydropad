  SUBROUTINE get_sp_emiss(rneh,Ti,sp_e)
!...returns spectrum of differential emissivity in photon/s/cm^3/keV
!
   USE xproj
   IMPLICIT NONE

   REAL*8, DIMENSION(:) :: sp_e(ibin1:ibin2)
   REAL*8 :: rneh,Ti
   INTEGER itemp1,itemp2
   REAL*8 T1,T2
!
  IF(Ti.ge.Tmin) THEN
!...particles hotter than the spectrum
    IF(Ti.ge.Tmax) THEN
      Ti=Tmax - 0.1
    ENDIF
!...setting suffixes of Power Coefficient files and related temperatures
    itemp1=(log10(Ti)-log10(Tmin))/dTT
    itemp2=itemp1+1
    T1=10**(log10(Tmin)+itemp1*dTT)
    T2=10**(log10(Tmin)+itemp2*dTT)
!
!...spectrum of differential emissivity in photon/s/cm^3/keV
    sp_e = 0.
    CALL get_pow(Ti,T1,T2,itemp1,itemp2,sp_e)
    sp_e = sp_e * rneh
!
  ENDIF
!
  END SUBROUTINE get_sp_emiss
!............................................................................
!
  SUBROUTINE get_pow(Ti,T1,T2,itemp1,itemp2,sp_e)
!...returns the differential spectrum in photon/s*cm^3/keV
!...of the power coefficient at the temperature Ti

  USE xproj
  IMPLICIT NONE

  INTEGER :: itemp1,itemp2
  REAL*8, DIMENSION(:) :: sp_e(ibin1:ibin2)
  REAL*8, ALLOCATABLE :: pow1(:),pow2(:)
  REAL*8 :: Ti,T1,T2
!
  ALLOCATE(pow1(ibin1:ibin2))
  ALLOCATE(pow2(ibin1:ibin2))
!
  pow1=pow(itemp1,ibin1:ibin2)
  pow2=pow(itemp2,ibin1:ibin2)
  WHERE(pow1.gt.0) 
    sp_e = exp( log(pow1) + (log(pow2)-log(pow1)) * &
                (log(Ti)-log(T1))/(log(T2)-log(T1)) )
  END WHERE
  WHERE(pow1.eq.0.and.pow2.gt.0) 
    sp_e = pow2*(log(Ti)-log(T1))/(log(T2)-log(T1))
  END WHERE
!
  DEALLOCATE(pow1)
  DEALLOCATE(pow2)
!
  END SUBROUTINE get_pow


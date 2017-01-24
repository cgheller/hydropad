MODULE HEAT
!
        IMPLICIT NONE
        SAVE
!
        REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: sncell
!
        REAL*8, PARAMETER :: epsilon_sn=3.1622777e25 ! ie: 10^25.5
        REAL*8, PARAMETER :: nu_sn=1.0/56.0
        REAL*8, PARAMETER :: fstar_sn=0.06
        REAL*8, PARAMETER :: fb_sn=0.08
!
! free parameters "eff", "eff1" (density threshold)
!
        REAL*8, PARAMETER :: eff=25.0
        REAL*8, PARAMETER :: eff1=5.0
!
END MODULE HEAT

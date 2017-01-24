MODULE xproj
!
	IMPLICIT NONE
	SAVE
!
! powd1 = number of considered temperature
! powd2 = number of frequency bins
!
	INTEGER, PARAMETER :: ntemp=91
	INTEGER, PARAMETER :: nchan=999
!
	REAL*8 :: powmin,powmax,de,Tmin,dTT,zmet,Tmax
	REAL*8 :: pow(0:ntemp-1,nchan)
!
	CHARACTER*90 :: filepow
!
	REAL*8 :: ez1,ez2
	INTEGER :: ibin1,ibin2
!
END MODULE xproj

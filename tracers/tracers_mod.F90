MODULE tracers_mod
!
	INTEGER ntrack
	REAL*8  zstarttrack
	INTEGER ntcluster
	INTEGER troutstep
	INTEGER checkstarttrack
	INTEGER size
	INTEGER winvx,winvy,winvz
	INTEGER first_time

	REAL*8, DIMENSION (:,:), ALLOCATABLE :: xtrack
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: ytrack
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: ztrack
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: xtrack0
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: ytrack0
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: ztrack0
	INTEGER, DIMENSION (:),  ALLOCATABLE :: xtstart
	INTEGER, DIMENSION (:),  ALLOCATABLE :: ytstart
	INTEGER, DIMENSION (:),  ALLOCATABLE :: ztstart
	INTEGER, DIMENSION (:),  ALLOCATABLE :: ntcell
	INTEGER, DIMENSION (:),  ALLOCATABLE :: ntrackers
        INTEGER, DIMENSION (:),  ALLOCATABLE :: njump

END MODULE tracers_mod

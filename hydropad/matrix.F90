MODULE matrix
!
USE dimension
IMPLICIT NONE
SAVE
!
! define baryonic and dark matter basic three dimensional vectors
!
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: nes3d,cho3d
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: p3d,rho3d,vx3d,vy3d,vz3d
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: p3dnew,rho3dnew,vx3dnew,vy3dnew,vz3dnew,cho3dnew
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: ttt,pold3d
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: vxold,vyold,vzold,rhoold
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: phi3d,phiold3d
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: rhodm3d
!
REAL*8, DIMENSION (:,:), ALLOCATABLE :: ppos
REAL*8, DIMENSION (:,:), ALLOCATABLE :: pvel
!
REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: grav_shap 
!
REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: gforce
!
!$acc declare create(p3d,rho3d,vx3d,vy3d,vz3d,gforce)
!
END MODULE MATRIX

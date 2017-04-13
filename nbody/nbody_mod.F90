MODULE nbody_mod
IMPLICIT NONE
SAVE

REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xsend,xrecv
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: vsend,vrecv
INTEGER, ALLOCATABLE, DIMENSION(:) :: nsendpe,nrecvpe,indexsend,psend,precv
INTEGER, ALLOCATABLE, DIMENSION(:) :: listofparticles
INTEGER :: nsend,nrecv

END MODULE nbody_mod

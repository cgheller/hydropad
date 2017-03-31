MODULE mpi_inc
!
SAVE
!
#ifdef USEMPI
INCLUDE 'mpif.h'
#endif
!
INTEGER :: ierr,npes
!
INTEGER :: ndims=3
INTEGER, ALLOCATABLE, DIMENSION(:) :: dims
LOGICAL, ALLOCATABLE, DIMENSION(:) :: periods
INTEGER, ALLOCATABLE, DIMENSION(:) :: coordinates
INTEGER :: reorder
INTEGER :: COMM_CART
INTEGER :: planexznbound, planeyznbound
! 1=x 1 left, 2 right; 2=y 3 front, 4 rear; 3=z 5 down, 6 up
INTEGER, ALLOCATABLE, DIMENSION(:) :: neighbour
INTEGER, ALLOCATABLE, DIMENSION(:) :: coords
#ifdef USEMPI
INTEGER :: status(MPI_STATUS_SIZE)
#endif
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: boxl
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: boxsize
INTEGER :: output_pe
!
END MODULE mpi_inc

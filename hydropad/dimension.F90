MODULE DIMENSION
!
! basic parameters to dimensionalize variables
!
IMPLICIT NONE
SAVE
!
! total mesh size
!
INTEGER :: ngridx, ngridy, ngridz
!
! mesh size per pe
!
INTEGER :: ngridxpe, ngridype, ngridzpe
!
! number of MPI ranks per dimension
!
INTEGER :: npesx, npesy, npesz
!
! number of tiles per dimension
!
INTEGER :: ntilex, ntiley, ntilez
!
! number of ghost cells per side per dimension
!
INTEGER :: nbound
!
! mesh size per pe with ghost regions
!
INTEGER :: nx,ny,nz
!
! total number of particles and number of particles per pe
!
INTEGER :: npart,npartpe
!
INTEGER :: ngrid
!
! OLD STUFF
!
INTEGER :: nxnynz
!
INTEGER :: n1
INTEGER :: n1n1,n11,n12,n121,n1pe,n11pe,nparmax
INTEGER :: ntot,n21,ngdim2,ngdim2pe,ngr2
INTEGER :: nxmax
!
INTEGER :: desca(12)
!
! rhotable -> numero di rette per il cooling dipendente da rho 
! rho2table -> numero di rette per il cooling dipendente da rho**2
! rho3table -> numero di rette per il cooling dipendente da rho**3
!
INTEGER :: rhotable, rho2table, rho3table
!
! number of bytes for each variable in output files
!
INTEGER, PARAMETER :: nlrec=4

END MODULE DIMENSION

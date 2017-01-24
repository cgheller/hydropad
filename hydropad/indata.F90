#include "hydrompi.def"
!       
SUBROUTINE indata(startdump)
!
	USE dimension
	USE matrix
	USE vector
	USE scalar
	USE io_mod
	USE tracers_mod
!
	USE mpi_inc
!
	IMPLICIT NONE
!
! define local variables
!
	INTEGER :: ii,iunit
        INTEGER :: i,j,k
	INTEGER :: startdump
	character*3 :: outform
	common/form/outform
!
! read simulation parameters
!
	  iunit=800+mype
	  OPEN(unit=iunit,file='input.dat',status='old')
!
! simulation related parameters
! explanation: 
!	       - cosmological_model : 60 characters string
!	       - owner : 60 characters string
!	       - date_of_production : ggmmaaaa
!	       - s_id : id number of the simulation
!	       - ngrid:	    number of cells per dimension (assuming cubic box)
!     	       - startdump: 0 starts from initial conditions
!                           1 starts from dump
!              - stopdump : 0 continue after dump
!                           1 stop after dump
!              - checkyn  : 0 do not check read data
!                           1 check read data
!              - ndump    : number of steps between two dumps
!              - tstop    : final redshift
!              - nuv      : use reheating (0=no heating
!                                          1=sn winds)
!              - ncool    : use cooling (1=true)
!	       - zstarttrack : tracers initial redshift
!              - box      : box size (Mpc/h)
!              - hnow     : Hubble parameter
!              - omega0   : density parameter
!              - dmfact   : DM density parameter
!              - bmfact   : BM density parameter
!              - lambda   : cosmological constant
!              - tempnow  : present BM average temperature
!              - hfrac    : hydrogen mass fraction
!	       - outform  : output data format
!              - tout     : output redshifts (max 9)
! 	    
!
	  read(iunit,'(a60)')cosmological_model
	  read(iunit,'(a60)')owner
	  read(iunit,*)date_of_production
	  read(iunit,*)s_id
	  read(iunit,*)ngrid
	  read(iunit,*)startdump
	  read(iunit,*)stopdump
	  read(iunit,*)checkyn
	  read(iunit,*)ndump
	  read(iunit,*)tstop_aux
	  read(iunit,*)nuv
	  read(iunit,*)ncool
	  read(iunit,*)zstarttrack
!
! physical parameters
!
	  read(iunit,*)box
	  read(iunit,*)hnow
	  read(iunit,*)dmfact
	  read(iunit,*)bmfact
	  read(iunit,*)lambda_vac
	  read(iunit,*)tempnow
	  read(iunit,*)hfrac
	  read(iunit,'(a3)')outform
!
! output redshifts (max 9)
!
	  do ii=1,9
	    read(iunit,*)tout(ii)
	  enddo
	  close(iunit)
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! define basic dimensional quantities
!
	npart=ngrid*ngrid*ngrid
        nx=ngrid
        ny=ngrid
        nz=ngrid/npes
!
	n1=ngrid
	nparmax=npart/npes
	nxmax=n1
        xmax=float(nxmax)
	n1n1=n1*n1
	n11=n1-1
	n12=n1/2
	n121=n12-1
	ntot = n1*n1*n1
	n21=n12+1
	ngr2=ntot/2
	ngdim2=ngr2+n1n1
        nxnynz=nx*ny*nz
        n1pe=n1/npes
	n11pe=n1pe-1

!
! allocate basic arrays
!
        CALL alloc_arrays
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! calculate fundamental density parameters:
!
	omega_dm=dmfact
	omega_bm=bmfact
	omega_m=omega_dm+omega_bm
	omega_lambda=lambda_vac
	omega_tot=omega_dm+omega_bm+omega_lambda
	omega_a=omega_m-0.3*omega_tot+0.3
	omega0=omega_tot
!
	cour=0.9
	ca=0.020
	eta1=0.3
!	eta2=0.005
	eta2=0.020
	dmax=2.0e10
	boltz=1.380658e-16
	mh=1.6605e-24
	norm=1.0
	gamma=1.66666667
!
! inversion of dmfact and bmfact
!
	if(dmfact.ne.0.0)then
	   dmfact=1.0/dmfact
	else
	   dmfact=1.0
	endif
	rdmfact=1.0/dmfact
	rbmfact=bmfact
	bmfact=1.0/bmfact
	box=box/hnow
	dxmpc=box/float(nx)
!
! calculate derivate quantities 
!
        amass=rdmfact*float(ntot)/float(npart)
	r=boltz/mh
	dx=1.0*norm
	rdx=1.0/dx
	rgamma=1.0/gamma
	gamma1=(gamma+1.0)/(gamma-1.0)
	rgamma1=1.0/(gamma-1)
	m=(gamma-1.0)/(2.0*gamma)
	rm=1.0/m
	gf=(gamma+1.0)/(2.0*gamma)
!
! initialize k space
!
        CALL start
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! initialize shape vector
!
	CALL grav_shape 
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!	
END SUBROUTINE indata

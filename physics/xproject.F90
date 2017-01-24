!
SUBROUTINE xproject
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
!
! local variables
!
	IMPLICIT NONE
	REAL*8, DIMENSION (:,:),ALLOCATABLE :: bright
	REAL*8, DIMENSION (:,:),ALLOCATABLE :: sz
	REAL*8, DIMENSION (:,:),ALLOCATABLE :: cndens
	REAL*8, DIMENSION (:,:),ALLOCATABLE :: ttt2d
	REAL*8, DIMENSION (:,:),ALLOCATABLE :: mat_aux
!
	CHARACTER*1  :: ch1,ch2,ch3,ch4,ch5
	CHARACTER*30 :: filename1
	INTEGER :: i1,j1,k1,nflag,nxny,nnzmax,kaux
	INTEGER :: i,j,k
	INTEGER :: ixmin,ixmax,iymin,iymax,izmin,izmax
	INTEGER :: i1aux,i2aux
	REAL*8  :: x1aux,x2aux,x3aux,x4aux,x5aux
	REAL*8  :: trad,sigmat,me,coh
	REAL*8  :: cfsz,cfemiss,cfcndens
	REAL*8  :: da,dl,dxf
	REAL*8  :: delta,cfflux,bcell,fbandflux
	REAL*8  :: fhy,fhe
	REAL*8  :: xnp,xne,rneh,tkev
	REAL*8  :: dxmpc_aux
!
! cfsz: electron mass (me) in kev
! Thomson cross section in cm**2. trad is radiation temperature in K
! The terms in (1+z) appears because rho is the comoving density
! 
	nxny=nx*ny
	fhy=hfrac
	fhe=1.0-hfrac
	if(redshiftold.eq.0.0)redshiftold=0.001
        trad=2.73
        sigmat=6.65e-25
        me=511.
	coh=2*2.998e3/hnow
!
! cfsz    = sz effect constant                 times
!           density dimensionalization const   times
!           Mpc to cm conversion factor        times
!           expansion contribute to comoving density
! 
        cfsz=sigmat*(fhy+0.5*fhe)*rhonow/mp*kb/me*(1.+redshiftold)**3
!
! cfemiss = free-free emission constant        times
!           density dimensionalization const   times
!           Mpc to cm conversion factor        times
!           expansion contribute to comoving density
!
	cfemiss=rhonow/mp*rhonow/mp*(1+redshiftold)**6
!
! cfcndens = column density dimensionalization factor
!
	cfcndens=(fhy+0.5*fhe)*rhonow/mp*(1.+redshiftold)**3
!
! compute diameter distance and luminosity distance: expression valid
! only for omega=1. Distances are espressed in Mpc.
!
	da=coh/(1.+redshiftold)*(1-1./sqrt(1.+redshiftold))
	dl=(1.+redshiftold)**2*da	
!
! 'dxmpc' is the grid spacing in Mpc comoving, 'dxf' is the corresponding 
! physical size in Mpc
!
	dxmpc_aux=dxmpc/(1.+redshiftold)
	dxf=dxmpc_aux*mpctocm
! 
! compute the angle which corresponds to the grid spacing. alpha gives the 
! same angle in arcsec
!  bright is in erg/sec/cm^2 (epsilon dx)
!  bright / 4 pi is the surface brightness from the cell erg/sec/cm^2/sterad
!  delta bright / 4 pi is the flux from the cell (delta is the angle 
!  under which the cell is seen)
!
	delta=dxmpc_aux/da
!	alpha=180.0*delta/pi*3600.
	alpha=(180.0*60/pi)
	cfflux=delta**2/4*pi
!
! compute surface brightness, sz effect and column density along 
! principal axis. The integration is performed over all the cells whose
! centre is within the virial radio of the corresponding cluster	  
!
!-------------------------------------------------------------------------
! integration along x-axis
!
	allocate(bright(ny,nz))
	allocate(sz(ny,nz))
	allocate(cndens(ny,nz))
	allocate(ttt2d(ny,nz))
	bright=0.0
	sz=0.0
	cndens=0.0
	ttt2d=0.0
!
	izmin=1
	izmax=nz
	iymin=1
	iymax=ny
!	     
	do k=izmin,izmax
	   k1=k
	   do j=1,ny
	      j1=j	     
	      do i=1,nx
	         i1=i
!		  		
!
		 xnp=fhy*rho3d(i1,j1,k1)
		 xne=(fhy+0.5*fhe)*rho3d(i1,j1,k1)
		 rneh=cfemiss*xnp*xne
!
		 tkev=ttt(i1,j1,k1)*kb
		 call flux(rneh,tkev,bcell)
!
                 bright(j1,k1)=bright(j1,k1)+bcell*dxf
                 sz(j1,k1)=sz(j1,k1)+rho3d(i1,j1,k1)*ttt(i1,j1,k1)*dxf
	         cndens(j1,k1)=cndens(j1,k1)+rho3d(i1,j1,k1)*dxf
	         ttt2d(j1,k1)=ttt2d(j1,k1)+ttt(i1,j1,k1)
	      enddo
! 		
              bright(j1,k1)=bright(j1,k1)*cfflux/(1.0+redshiftold)**4
! SZ effect
	      sz(j1,k1)=sz(j1,k1)*cfsz
! optical depth
	      cndens(j1,k1)=cfcndens*cndens(j1,k1)
! average volume weighted temperature
	      ttt2d(j1,k1)=ttt2d(j1,k1)/float(nx)
!	   
           enddo
        enddo 
!
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(mype.eq.0)then
	  write(*,*) 'End of integration along x direction'
	endif
!
101     format(2(1x,i4),5(1x,e11.5))
!
	call getindex10000(nstep,ch1,ch2,ch3,ch4,ch5)
!
	filename1='proj_x.'//ch1//ch2//ch3//ch4//ch5//'.dat'
	open(201,file=filename1)
	do i=0,npes-1
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(mype.eq.i)then
	nnzmax=ny*nz*mype
	do k=1,nnzmax
	  read(201,*)i1aux,i2aux,x1aux,x2aux,x3aux,x4aux,x5aux
	enddo
	do k=1,nz
	kaux=k+mype*nz
	do j=1,ny

! write order : flux, SZ, column density, avg temperature, surface brightness

	  write(201,101)j,kaux,bright(j,k),sz(j,k),cndens(j,k),&
                        ttt2d(j,k),bright(j,k)/delta**2/alpha**2
	enddo	
	enddo	
	close(201)
	endif
	enddo
!
	deallocate(bright)
	deallocate(sz)
	deallocate(cndens)
	deallocate(ttt2d)
!      
!-------------------------------------------------------------------------
! integration along y-axis
!
	allocate(bright(nx,nz))
	allocate(sz(nx,nz))
	allocate(cndens(nx,nz))
	allocate(ttt2d(nx,nz))
	bright=0.0
	sz=0.0
	cndens=0.0
	ttt2d=0.0
!
	do k=izmin,izmax
	   k1=k
	   do i=1,nx
	      i1=i
	      do j=1,ny
	         j1=j
!
!
		 xnp=fhy*rho3d(i1,j1,k1)
		 xne=(fhy+0.5*fhe)*rho3d(i1,j1,k1)
		 rneh=cfemiss*xnp*xne
!
		 tkev=ttt(i1,j1,k1)*kb
		 call flux(rneh,tkev,bcell)
!
                 bright(i1,k1)=bright(i1,k1)+bcell*dxf
                 sz(i1,k1)=sz(i1,k1)+rho3d(i1,j1,k1)*ttt(i1,j1,k1)*dxf
	         cndens(i1,k1)=cndens(i1,k1)+rho3d(i1,j1,k1)*dxf
	         ttt2d(i1,k1)=ttt2d(i1,k1)+ttt(i1,j1,k1)
  	      enddo
! 		
! Surface brightness
!
              bright(i1,k1)=bright(i1,k1)*cfflux/(1.0+redshiftold)**4
! SZ effect
	      sz(i1,k1)=sz(i1,k1)*cfsz
! optical depth
	      cndens(i1,k1)=cfcndens*cndens(i1,k1)
! average volume weighted temperature
	      ttt2d(i1,k1)=ttt2d(i1,k1)/float(ny)
           enddo
	enddo 
!
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(mype.eq.0)then
	  write(*,*) 'End of integration along y direction'
	endif
!
	filename1='proj_y.'//ch1//ch2//ch3//'.dat'
	open(201,file=filename1)
	do i=0,npes-1
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(mype.eq.i)then
	nnzmax=ny*nz*mype
	do k=1,nnzmax
	  read(201,*)i1aux,i2aux,x1aux,x2aux,x3aux,x4aux,x5aux
	enddo
	do k=1,nz
	kaux=k+mype*nz
	do j=1,ny
	  write(201,101)j,kaux,bright(j,k),sz(j,k),cndens(j,k),&
                        ttt2d(j,k),bright(j,k)/delta**2/alpha**2
	enddo
	enddo
        close(201)
	endif
	enddo
!
	deallocate(bright)
	deallocate(sz)
	deallocate(cndens)
	deallocate(ttt2d)
!      
!-------------------------------------------------------------------------
! integration along z-axis
!
	allocate(bright(nx,ny))
	allocate(sz(nx,ny))
	allocate(cndens(nx,ny))
	allocate(ttt2d(nx,ny))
	allocate(mat_aux(nx,ny))
!
	bright=0.0
	sz=0.0
	cndens=0.0
	ttt2d=0.0
	mat_aux=0.0
!
        do j=1,ny
	   j1=j
           do i=1,nx
	      i1=i
!
	      do k=1,nz
	         k1=k
!
!
		 xnp=fhy*rho3d(i1,j1,k1)
		 xne=(fhy+0.5*fhe)*rho3d(i1,j1,k1)
		 rneh=cfemiss*xnp*xne
!
		 tkev=ttt(i1,j1,k1)*kb
		 call flux(rneh,tkev,bcell)
!
                 bright(i1,j1)=bright(i1,j1)+bcell*dxf
                 sz(i1,j1)=sz(i1,j1)+rho3d(i1,j1,k1)*ttt(i1,j1,k1)*dxf
	         cndens(i1,j1)=cndens(i1,j1)+rho3d(i1,j1,k1)*dxf
	         ttt2d(i1,j1)=ttt2d(i1,j1)+ttt(i1,j1,k1)
	      enddo
! 		
! Surface brightness
!
              bright(i1,j1)=bright(i1,j1)*cfflux/(1.0+redshiftold)**4
! SZ effect
	      sz(i1,j1)=sz(i1,j1)*cfsz
! optical depth
	      cndens(i1,j1)=cfcndens*cndens(i1,j1)
! average volume weighted temperature
	      ttt2d(i1,j1)=ttt2d(i1,j1)/float(nz)
           enddo
	enddo 
!
	CALL MPI_Reduce(bright,mat_aux,nxny,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
	bright=mat_aux
	mat_aux=0.0
	CALL MPI_Reduce(sz,mat_aux,nxny,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
	sz=mat_aux
	mat_aux=0.0
	CALL MPI_Reduce(cndens,mat_aux,nxny,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
	cndens=mat_aux
	mat_aux=0.0
	CALL MPI_Reduce(ttt2d,mat_aux,nxny,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
	ttt2d=mat_aux/float(npes)

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if(mype.eq.0)then
	  write(*,*) 'End of integration along z direction'
	endif
!
	if(mype.eq.0)then
	filename1='proj_z.'//ch1//ch2//ch3//'.dat'
	open(201,file=filename1)
        do k=1,ny
        do j=1,nx
	  write(201,101)j,k,bright(j,k),sz(j,k),cndens(j,k),&
                        ttt2d(j,k),bright(j,k)/delta**2/alpha**2
        enddo
        enddo
        close(201)
	endif
!
	deallocate(bright)
	deallocate(sz)
	deallocate(cndens)
	deallocate(ttt2d)
	deallocate(mat_aux)
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
END SUBROUTINE xproject

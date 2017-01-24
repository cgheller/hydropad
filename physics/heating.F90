
SUBROUTINE heating
!
! SUPERNOVAE SHOCK HEATING
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
!
	USE heat
!
! local variables
!
	IMPLICIT NONE
	integer :: req(4),status_array(MPI_STATUS_SIZE,4)
	integer :: nslice
	integer :: i,j,k
	integer :: from_pe,to_pe
	integer :: indp1,indm1
	integer :: ii,jj,kk,i_index,j_index,k_index,nrst
	integer :: imax,jmax,kmax,imin,jmin,kmin
	integer :: nnv,nnjeans,nndyn
	integer :: nnv0,nnjeans0,nndyn0
	real*8  :: vup,vdown
	real*8  ::somma,mb
	real*8  :: gradvx,gradvy,gradvz,div_v
	real*8  :: rhon,rhoe,vsound1,deltabm,deltadm
	real*8  :: rhotot,mjfact,mjeans,tdyn,rtcool,tcool
	real*8  :: fh,lambda,lambda_rad,rtdyn
        real*8  :: lambda_com,lambda_tot,eaux
	real*8  :: tsn,tsn_d,e_inputz,rst,sigmarst
	real*8  :: dist,av_e,e_weight,xcen,ycen,zcen
	real*8  :: mcell,mstar
	real*8  :: tttaux
	real*8  :: fact_sn,fact_stars
	real*8  :: xn_stars,xn_sn,xnn_sn,xnn_sn_tot,xnn_sn_cell
	real*8  :: snmass,snmasspe,metal
	real*8  :: time8,timew
	real*8  :: star_formation_rate,star_formation_rate_tot
        real*8, dimension (:,:,:), allocatable :: vz_aux
	real*8, parameter :: eta_sn=5.5
!
        allocate (vz_aux(nx,ny,2))
	vz_aux=0.0
	nslice=nx*ny
	if(nstep.eq.1)ntotsn=0.0
!
        from_pe=mype-1
        if(from_pe.lt.0)from_pe=npes-1
        to_pe=mype+1
        if(to_pe.gt.npes-1)to_pe=0
        CALL MPI_Irecv(vz_aux(1,1,1),nslice,MPI_DOUBLE_PRECISION,&
                       from_pe,30,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_Isend(vz3d(1,1,nz),nslice,MPI_DOUBLE_PRECISION,&
                       to_pe,30,MPI_COMM_WORLD,req(2),ierr)
        CALL MPI_Irecv(vz_aux(1,1,2),nslice,MPI_DOUBLE_PRECISION,&
                       to_pe,40,MPI_COMM_WORLD,req(3),ierr)
        CALL MPI_Isend(vz3d(1,1,1),nslice,MPI_DOUBLE_PRECISION,&
                       from_pe,40,MPI_COMM_WORLD,req(4),ierr)
!
        CALL MPI_WAITALL(4,req,status_array,ierr)
!
	time8=3.1536e15/tnow
	timew=dt/time8
!
! calculate energy release per unit comoving volume 
!
! e = a^3 * E / V
!
! e = a^3 * b * M / V = b * ( a^3 * M / V ) = b * ( rho_code * rho0 )
!
! e_code * rho0 * x0^2 / t0^2 = b * rho_code * rho0
!
! e_code = b * rho_code * t0^2 / x0^2
!
! b = eff * epsilon_sn * nu_sn * fstar_sn * fb_sn / M_sun * dt8
!
! dt8 and M_sun are required to get the correct units
! Notice that the mass in stars = mass in supernovae / nu_sn
!
        fact_sn=eff*tsux*epsilon_sn*nu_sn*fstar_sn*fb_sn*&
                (epsilon_sn/msun)
!
        fact_sn=timew*fact_sn
!
! calculate the mass in supernovae
!
! Msn = eps * c * M
!
! Msn = eps * c * rho0 * rho_code / a^3 * dx_code^3 * x0^3 * a^3
!
! Msn = eps * c * rho0 * rho_code * x0^3
!
! eps = dt8, c = eff * nu_sn * fstar_sn * fb_sn
!
	xn_stars=eff*fstar_sn*fb_sn*(rhonow*xnow*xnow)*&
              (xnow/msun)
        xn_sn=nu_sn*xn_stars*timew
!
        fh=0.75
!
	nnv=0
	nnjeans=0
	nndyn=0

        do k=1,nz
        do j=1,ny
        do i=1,nx
	  nes3d(i,j,k)=0
!
! check overdensity
!
	  if(rho3d(i,j,k).ge.eff1*rbmfact)then
!
! calculate velocity divergence
!
	  indp1=i+1
	  if(indp1.gt.nx)indp1=1
	  indm1=i-1
	  if(indm1.lt.1)indm1=nx
	  gradvx=vx3d(indp1,j,k)-vx3d(indm1,j,k)
!
	  indp1=j+1
	  if(indp1.gt.ny)indp1=1
	  indm1=j-1
	  if(indm1.lt.1)indm1=ny
	  gradvy=vy3d(i,indp1,k)-vy3d(i,indm1,k)
!
	  indp1=k+1
	  if(indp1.le.nz)vup=vz3d(i,j,indp1)
	  if(indp1.gt.nz)vup=vz_aux(i,j,2)
	  indm1=k-1
	  if(indm1.ge.1)vdown=vz3d(i,j,indm1)
	  if(indm1.lt.1)vdown=vz_aux(i,j,1)
	  gradvz=vup-vdown
!
	  div_v=gradvx+gradvy+gradvz
	  if(div_v.gt.0.0)then
	    goto 333
	  else
	    nnv=nnv+1
	  endif
!
! calculate jeans mass and baryonic mass of a cell
!
	  rhotot=rho3d(i,j,k)+rhodm3d(i,j,k)
	  mjeans=(pi/(gconst*rhotot))**1.5
          vsound1=sqrt(gamma*p3d(i,j,k)/rho3d(i,j,k))
	  deltabm=rho3d(i,j,k)-1.0
	  deltadm=rhodm3d(i,j,k)-1.0
	  mjeans=mjeans*vsound1**3*rho3d(i,j,k)
	  mb=rho3d(i,j,k)
          if(mb.lt.mjeans)then
	    goto 333
	  else
	    nnjeans=nnjeans+1
	  endif
!
! calculate dynamical time and cooling time of a cell
!
	  rhon=rho3d(i,j,k)
	  tdyn=sqrt(3.0*pi/(32.0*gconst*rhotot))
!
	  tttaux=ttt(i,j,k)
	  lambda_rad=hfact*rhon**2*lambda(tttaux)/at**3
	  lambda_com=cmfact*5.41d-36*(1.0+redshift)**4*&
                     ttt(i,j,k)*rhoe(rhon,ttt(i,j,k),fh)
	  lambda_tot=lambda_rad+lambda_com
!
	  eaux=p3d(i,j,k)/(gamma-1.0)
	  rtcool=lambda_tot/eaux
	  rtdyn=1.0/tdyn
          if(rtcool.le.rtdyn)then
            goto 333
	  else
	    nndyn=nndyn+1
	  endif
!
! if star formation conditions are true nes3d=1
!
	  nes3d(i,j,k)=1
!
333	  continue
	  endif
	enddo
	enddo
	enddo
!
	nnv0=0
	nnjeans0=0
	nndyn0=0
!
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
	CALL MPI_Allreduce(nnv,nnv0,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
	CALL MPI_Allreduce(nnjeans,nnjeans0,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
	CALL MPI_Allreduce(nndyn,nndyn0,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
!
	if(mype.eq.0)&
           write(500,89)redshift,nnv0,nnjeans0,nndyn0
89	format(1x,e11.5,3(1x,i8))
!
	xnn_sn=0.0
	star_formation_rate=0.0
!
        do k=1,nz
        do j=1,ny
        do i=1,nx
	  e_inputz=0.0
	  somma=0.0
	  if(nes3d(i,j,k).eq.1)then
!
! calculate SN TOTAL energy release of a cell (in code units)
!
	    rhotot=rho3d(i,j,k)+rhodm3d(i,j,k)
            mcell=rho3d(i,j,k)
!
! energy input (code units)
!
            e_inputz=fact_sn*mcell
!
! xnn_sn is the fraction of baryonic mass of the cell that ends in supernovae
!
! mass of the cell in supernovae (solar masses)
!
	    xnn_sn_cell=xn_sn*mcell
!
! time unit for the star formation rate: 10^8 years
!
	    star_formation_rate=star_formation_rate+xn_stars*mcell
!
! total mass in supernovae (solar masses)
!
            xnn_sn=xnn_sn_cell+xnn_sn
	    sncell(i,j,k)=sncell(i,j,k)+xnn_sn_cell
!
! calculate Sedov-Taylor radius (in code units)
!
	    rst=0.76*(e_inputz*dt**2/mcell)**0.2/at
	    if(rst.gt.1.0)then
	       write(50,*)redshift,rst
	    endif
!
! at the moment Sedov-Taylor spatial distribution is NOT implemented.
! in fact we set:
!
	    rst=1.0
!
	    if(rst.le.1.0)then
	      e_weight=e_inputz*(gamma-1.0) 
	      p3d(i,j,k)=p3d(i,j,k)+e_weight
	    else
	      xcen=float(i)-0.5
	      ycen=float(j)-0.5
	      zcen=float(k)-0.5
!
! determine cells inside box of radius 2*rst
!
	      nrst=nint(rst)
	      imin=i-nrst
	      jmin=j-nrst
	      kmin=k-nrst
	      imax=i+nrst
	      jmax=i+nrst
	      kmax=k+nrst
	      somma=0.0
	      do kk=kmin,kmax
	      k_index=kmin
	      if(k_index.le.0)k_index=k_index+ngrid
	      if(k_index.gt.ngrid)k_index=k_index-ngrid
	      do jj=jmin,jmax
	      j_index=jmin
	      if(j_index.le.0)j_index=j_index+ngrid
	      if(j_index.gt.ngrid)j_index=j_index-ngrid
	      do ii=imin,imax
	      i_index=imin
	      if(i_index.le.0)i_index=i_index+ngrid
	      if(i_index.gt.ngrid)i_index=i_index-ngrid
!
	      dist=sqrt((xcen-(ii-0.5))**2+&
                        (ycen-(jj-0.5))**2+&
                        (zcen-(kk-0.5))**2)
!
! calculate pressure correction term
!
	      if(dist.lt.rst)then
	        sigmarst=rst/3.0
	        av_e=(gamma-1.0)*e_inputz/(sqrt(pi)*sigmarst)**3
	        e_weight=av_e*exp(-dist**2/sigmarst**2)
	        somma=somma+e_weight
!
!!HPF$ CRITICAL
!
!	        p3d(i_index,j_index,k_index)=
!                p3d(i_index,j_index,k_index)+e_weight
!!HPF$ END CRITICAL
	      endif
!
	      enddo
	      enddo
	      enddo
	    endif
!
! end of the Sedorf-Taylor distribution part (NOT YET IMPLEMENTED)
!
! update temperature
!
            ttt(i,j,k)=tfact*p3d(i,j,k)/rho3d(i,j,k)
	  endif
112     continue
	enddo
	enddo
	enddo
!
! calculate total number of supernovae produced since the beginning
!
	xnn_sn_tot=0.0
	star_formation_rate_tot=0.0
        CALL MPI_Allreduce(xnn_sn,xnn_sn_tot,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
	ntotsn=ntotsn+xnn_sn_tot
        CALL MPI_Allreduce(star_formation_rate,star_formation_rate_tot,&
                           1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
	write(501,100)redshift,ntotsn,star_formation_rate_tot
100	format(3(1x,e13.7))
!
	deallocate (vz_aux)
!
END SUBROUTINE heating
!
!*********************************************************************
!*********************************************************************
!
FUNCTION lambda(t)
!
	IMPLICIT NONE
	real*8 :: lambda,t,tlog,exp1,exp2
!
	tlog=log10(t)
!
	exp1=-0.1-1.88*(5.23-tlog)**4
	if(tlog.le.6.2)then
	  exp2=-1.7-0.2*(6.2-tlog)**4
	else
	  exp2=-1.7
	endif
!
	lambda=10.0**(-21)*(10.0**exp1+10.0**exp2)
!
END FUNCTION lambda
!
!*********************************************************************
!*********************************************************************
!
FUNCTION rhoe(rho,t,fh)
!
	IMPLICIT NONE
	REAL*8 :: rho,t,rhoe,fh,ff
!
	ff=0.5*(fh+1)
!
	if(t.lt.1e4.and.t.ge.1.0)then
	  rhoe=ff*rho*exp(1.0-1e4/t)
        else if(t.lt.1)then
          rhoe=0.0
	else
	  rhoe=ff*rho
	endif
!
END FUNCTION rhoe

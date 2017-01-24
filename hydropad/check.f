	subroutine check
C
        include 'dim.h'
        include 'matrix.h'
        include 'vector.h'
        include 'scalar.h'
	include 'mpi_inc.h'
C
C local variables
C
	integer k_pe,iunit
C
	iunit=1400+mype
	k_pe=1
	write(iunit,*)nes3d(1,1,k_pe)
	write(iunit,*)cho3d(1,1,k_pe)
	write(iunit,*)p3d(1,1,k_pe)
	write(iunit,*)rho3d(1,1,k_pe)
	write(iunit,*)vx3d(1,1,k_pe)
	write(iunit,*)vy3d(1,1,k_pe)
	write(iunit,*)vz3d(1,1,k_pe)
	write(iunit,*)phi3d(1,1,k_pe)
	write(iunit,*)phiold3d(1,1,k_pe)
	write(iunit,*)ttt(1,1,k_pe)
	write(iunit,*)rhodm3d(1,1,k_pe)
	write(iunit,*)pold3d(1,1,k_pe)
	write(iunit,*)vxold(1,1,k_pe)
	write(iunit,*)vyold(1,1,k_pe)
	write(iunit,*)vzold(1,1,k_pe)
	write(iunit,*)rhoold(1,1,k_pe)
	write(iunit,*)x1(k_pe)
	write(iunit,*)x2(k_pe)
	write(iunit,*)x3(k_pe)
	write(iunit,*)v1(k_pe)
	write(iunit,*)v2(k_pe)
	write(iunit,*)v3(k_pe)
	write(iunit,*)grav_shap(0,0,0)
C
	write(iunit,*)nes(1)
	write(iunit,*)cho(1)
	write(iunit,*)pres(1)
	write(iunit,*)rho(1)
	write(iunit,*)vx(1)
	write(iunit,*)vy(1)
	write(iunit,*)vz(1)
	write(iunit,*)ghalf(1)
	write(iunit,*)g(1)
	write(iunit,*)c(1)
	write(iunit,*)eint(1)
	write(iunit,*)etot(1)
	write(iunit,*)rhoold1d(1)
	write(iunit,*)vxold1d(1)
	write(iunit,*)phi1(1)
	write(iunit,*)phi0(1)
C
	write(iunit,*)nstep,npl
	write(iunit,*)dt,t,told,at,dat,ath,dtold,dath
	write(iunit,*)rdtath,redshift,redshiftold,tin,tstop
	write(iunit,*)tstop_aux,atnew,thalf,rat,ca
	write(iunit,*)tout(1)
	write(iunit,*)bmfact,dmfact,rbmfact,rdmfact,dmax,amass
	write(iunit,*)gamma,gamma1,m,gf,re
	write(iunit,*)rm,rgamma,rgamma1
	write(iunit,*)cour,vsound,velmax
	write(iunit,*)eta1,eta2
	write(iunit,*)dx,rdx,norm,xmax
	write(iunit,*)hnow,xnow,tnow,rhonow,tempnow
	write(iunit,*)pfact,vnor,tfact,mfact,tsux,cmfact
	write(iunit,*)gconst,boltz,r,light
	write(iunit,*)mh,hfract,hfact,hfrac
	write(iunit,*)omega,omega0
	write(iunit,*)box
	write(iunit,*)ifail
	write(iunit,*)h_sec
	write(iunit,*)ak(0),akq(0)
	write(iunit,*)ndump,nuv
	write(iunit,*)stopdump
	write(iunit,*)pi,twopi
	write(iunit,*)mype
C
	return
	end

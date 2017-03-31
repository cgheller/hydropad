MODULE scalar
!
!USE dimension
IMPLICIT NONE
SAVE
!
! define basic scalar variables
!
! time related varibles:
!
CHARACTER(len=30)variablefilename
INTEGER :: maxsteps
INTEGER :: nstep,npl
INTEGER :: ifail
INTEGER :: startdump
INTEGER :: chooseinit
REAL(KIND=8) :: dt,t,told,at,dat,ath,dtold,dath,dth,dtinit
REAL(KIND=8) :: rdtath,redshift,redshiftold,tin,tstop,dtat
REAL(KIND=8) :: tstop_aux,atnew,datnew,thalf,rat,ca
REAL(KIND=8) :: tout(9)
REAL(KIND=8) :: xsinh
REAL(KIND=8) :: t0h0
REAL(KIND=8) :: h_sec,tnow_norm
REAL(KIND=8) :: hubble
REAL(KIND=8) :: zstart
!
! density related variables:
!
REAL(KIND=8) :: lambda_vac,bmfact,dmfact,rbmfact,rdmfact,dmax,amass
!
! fluid dynamics related variables:
!
REAL(KIND=8) :: gamma,gamma1,m,gf,re
REAL(KIND=8) :: rm,rgamma,rgamma1
REAL(KIND=8) :: cour,vsound,velmax
REAL(KIND=8) :: eta1,eta2
!
! grid related variables:
!
REAL(KIND=8) :: dx,rdx,norm,xmax
!
!
! normalization related variables
!
REAL(KIND=8) :: hnow,xnow,tnow,rhonow,tempnow
REAL(KIND=8) :: pfact,vnor,tfact,mfact,tsux,cmfact
!
! physical variables:
!
REAL(KIND=8) :: gconst,boltz,r,light
REAL(KIND=8) :: mh,hfract,hfact,hfrac
!
! cosmological variables:
!
REAL(KIND=8) :: omega,alpha,omega0
REAL(KIND=8) :: box,boxMpc_over_h
REAL(KIND=8) :: omega_dm,omega_bm,omega_m,omega_lambda
REAL(KIND=8) :: omega_tot,omega_a
!
! kspace vectors
!
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ak,akq
!
! input data
!
INTEGER :: ndump,nuv,ncool
INTEGER :: stopdump,checkyn
!
! constants
!
INTEGER :: mype,my_pe
!
REAL(KIND=8), parameter :: pi=3.1415926, twopi=2.0*pi
!
INTEGER :: noutput
!
REAL(KIND=8) :: year
REAL(KIND=8), parameter :: mpctocm=3.086e24
REAL(KIND=8), parameter :: year_in_secs=1.0e-6/31536000.0
REAL(KIND=8), parameter :: msun=1.9889e33
!
REAL(KIND=8) :: dxmpc
!
INTEGER :: nrot
!
REAL(KIND=8), PARAMETER :: mp=1.672623e-24
REAL(KIND=8), PARAMETER :: gaunt=1.2
REAL(KIND=8), PARAMETER :: kb=8.617e-8
!
REAL(KIND=8) :: ntotsn
!
REAL(KIND=8) :: zpreheating, temppreheating
REAL(KIND=8) :: enzoactive
!
! flattening parameters
!
REAL(KIND=8) :: mingradflat,flatvalue
!$acc declare copyin(mingradflat,flatvalue,gamma,rat,m,rgamma1,rm,rgamma,dmax,rdx,rdtath,dt,dx,dat,at,&
!$acc &              gamma1,gf,eta1,eta2)
END MODULE scalar

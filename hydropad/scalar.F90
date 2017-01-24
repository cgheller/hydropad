MODULE scalar
!
USE dimension
IMPLICIT NONE
SAVE
!
! define basic scalar variables
!
! time related varibles:
!
CHARACTER(len=30)variablefilename
integer :: maxsteps
integer :: nstep,npl
integer :: ifail
real*8 :: dt,t,told,at,dat,ath,dtold,dath,dth
real*8 :: rdtath,redshift,redshiftold,tin,tstop
real*8 :: tstop_aux,atnew,datnew,thalf,rat,ca
real*8 :: tout(9)
real*8 :: xsinh
real*8 :: t0h0
real*8 :: h_sec,tnow_norm
real*8 :: hubble
real*8 :: zstart
!
! density related variables:
!
real*8 :: lambda_vac,bmfact,dmfact,rbmfact,rdmfact,dmax,amass
!
! fluid dynamics related variables:
!
real*8 :: gamma,gamma1,m,gf,re
real*8 :: rm,rgamma,rgamma1
real*8 :: cour,vsound,velmax
real*8 :: eta1,eta2
!
! grid related variables:
!
real*8 :: dx,rdx,norm,xmax
!
!
! normalization related variables
!
real*8 :: hnow,xnow,tnow,rhonow,tempnow
real*8 :: pfact,vnor,tfact,mfact,tsux,cmfact
!
! physical variables:
!
real*8 :: gconst,boltz,r,light
real*8 :: mh,hfract,hfact,hfrac
!
! cosmological variables:
!
real*8 :: omega,alpha,omega0
real*8 :: box
real*8 :: omega_dm,omega_bm,omega_m,omega_lambda
real*8 :: omega_tot,omega_a
!
! kspace vectors
!
real*8, DIMENSION(:), ALLOCATABLE :: ak,akq
!
! input data
!
integer :: ndump,nuv,ncool
integer :: stopdump,checkyn
!
! constants
!
integer :: mype,my_pe
!
real*8, parameter :: pi=3.1415926, twopi=2.0*pi
!
integer :: noutput
!
real*8 :: year
real*8, parameter :: mpctocm=3.086e24
real*8, parameter :: year_in_secs=1.0e-6/31536000.0
real*8, parameter :: msun=1.9889e33
!
real*8 :: dxmpc
!
integer :: nrot
!
REAL*8, PARAMETER :: mp=1.672623e-24
REAL*8, PARAMETER :: gaunt=1.2
REAL*8, PARAMETER :: kb=8.617e-8
!
REAL*8 :: ntotsn
!
real*8 :: zpreheating, temppreheating
real*8 :: enzoactive
!
END MODULE scalar

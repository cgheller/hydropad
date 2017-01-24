#include "hydrompi.def"
!	
! NB: in this first impementation the 1D size MUST be divisible by the number 
! of MPI ranks in that dimension!!!
!
PROGRAM hydropad
!
USE dimension
USE matrix
USE vector
USE scalar
USE mpi_inc
!USE io_mod
#ifdef PROFILE
USE profile
#endif
IMPLICIT NONE
!
! local variables
!
INTEGER :: i,j,k
!
#ifdef PROFILE
t_cpu=0.0
#endif
!
! INITIALIZATION - MPI setup
!
#ifndef USEMPI
npes = 1
mype = 0
npesx = 1
npesy = 1
npesz = 1
write(*,*)
write(*,*)'The code is compiled without MPI directives'
write(*,*)
#endif
#ifdef USEMPI
!
! Create Default Communicator
!
CALL MPI_Init(ierr)
CALL MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
CALL MPI_Comm_rank(MPI_COMM_WORLD,mype,ierr)
#endif
!
! Initialize the simulation
!
CALL init_sys
!
! THIS HAS TO BE ALL REMOVED
!

if(mype .EQ. 0) then
write(*,*)ngridx
write(*,*)ngridy
write(*,*)ngridz
write(*,*)ngridxpe
write(*,*)ngridype
write(*,*)ngridzpe
write(*,*)npesx
write(*,*)npesy
write(*,*)npesz
write(*,*)nbound
endif
! END REMOVE
!

#ifdef PROFILE
if(mype.eq.0)then
  write(*,*)
  write(*,*)'WARNING ** !! ** the code is compiled '
  write(*,*)'for profiling tests'
  write(*,*)
endif
#endif
!
if(mype.eq.0)write(*,*)'*** STARTING SIMULATION ***'

!
! local variables -->
!
nstep=0
npl=1
!
! TIME LOOP STARTS
!
!!!!!!!do while(told.lt.tstop)
do while(nstep.lt.maxsteps)
!
  nstep=nstep+1
!
  if(mype.eq.0)then
    year=told*tnow*year_in_secs
    write(*,*)'Calculating step : ',nstep
    write(*,99)'t_i, t_f, dt : ',told,t,dt
99  format(1x,a31,3(1x,e13.7))
  endif
!
! system evolution
!
  CALL evolve_sys
  write(*,*)maxval(rho3d),minval(rho3d)
!
#ifdef USEMPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! TIME LOOP ENDS
!
enddo
!
! CLA: Do a very stupid I/O
#ifdef SHOCKX
k=nz/2
j=ny/2
do i = nbound+1,nx-nbound
  write(500,'(i3,2x,5(e13.7,1x))')i-nbound,rho3d(i,j,k),p3d(i,j,k),vx3d(i,j,k),vy3d(i,j,k),vz3d(i,j,k)
enddo
#endif
#ifdef SHOCKY
i=nx/2
k=nz/2
!j=ny/2
do j = nbound+1,ny-nbound
  write(500,'(i3,2x,5(e13.7,1x))')j-nbound,rho3d(i,j,k),p3d(i,j,k),vx3d(i,j,k),vy3d(i,j,k),vz3d(i,j,k)
enddo
#endif
#ifdef SHOCKZ
i=nx/2
j=ny/2
do k = nbound+1,nz-nbound
  write(500,'(i3,2x,5(e13.7,1x))')k-nbound,rho3d(i,j,k),p3d(i,j,k),vx3d(i,j,k),vy3d(i,j,k),vz3d(i,j,k)
enddo
#endif
#ifdef SEDOV
k=nz/2
j=ny/2
do i = nbound+1,nx-nbound
  write(500,'(i3,2x,5(e13.7,1x))')i-nbound,rho3d(i,j,k),p3d(i,j,k),vx3d(i,j,k),vy3d(i,j,k),vz3d(i,j,k)
enddo
#endif

variablefilename="rho.bin"
call write_var_mpi(rho3d,nx,ny,nz,variablefilename,mype)
variablefilename="pressure.bin"
call write_var_mpi(p3d,nx,ny,nz,variablefilename,mype)
variablefilename="vx.bin"
call write_var_mpi(vx3d,nx,ny,nz,variablefilename,mype)
variablefilename="vy.bin"
call write_var_mpi(vy3d,nx,ny,nz,variablefilename,mype)
variablefilename="vz.bin"
call write_var_mpi(vz3d,nx,ny,nz,variablefilename,mype)

!
! deallocate arrays
!
CALL dealloc_arrays

#ifdef USEMPI
CALL MPI_FINALIZE(ierr)
#endif
!
stop
END PROGRAM hydropad

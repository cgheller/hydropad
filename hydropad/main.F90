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
#ifdef PGIACC
USE cudafor
#endif
!USE io_mod
#ifdef PROFILE
USE profile
#endif
IMPLICIT NONE
!
! local variables
!
REAL(KIND=8) :: tstart,tend
INTEGER :: i,j,k,istat
#ifdef PGIACC
INTEGER(kind=cuda_count_kind) :: hsize
!
hsize = 32*1024*1024
!istat = cudaDeviceSetLimit(cudaLimitMallocHeapSize, hsize)
#endif
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
if(mype == 0)write(*,*)
if(mype == 0)write(*,*)'The code is compiled without MPI directives'
if(mype == 0)write(*,*)
#endif
#ifdef USEMPI
!
! Create Default Communicator
!
CALL MPI_Init(ierr)
CALL MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
CALL MPI_Comm_rank(MPI_COMM_WORLD,mype,ierr)
tstart = MPI_Wtime()
#endif
!
! Initialize the simulation
!
!$acc data present(rho3d,p3d,vx3d,vy3d,vz3d,gforce)
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
do while(nstep.lt.maxsteps .AND. told.lt.tstop)
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
!
#ifdef USEMPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! TIME LOOP ENDS
!
enddo
!
!$acc end data

#ifdef HYDRO
variablefilename="rho.bin"
CALL write_var_mpi(rho3d,nx,ny,nz,variablefilename,mype)
variablefilename="pressure.bin"
CALL write_var_mpi(p3d,nx,ny,nz,variablefilename,mype)
variablefilename="vx.bin"
CALL write_var_mpi(vx3d,nx,ny,nz,variablefilename,mype)
variablefilename="vy.bin"
CALL write_var_mpi(vy3d,nx,ny,nz,variablefilename,mype)
variablefilename="vz.bin"
CALL write_var_mpi(vz3d,nx,ny,nz,variablefilename,mype)
#endif
#ifdef NBODY
variablefilename="positions.bin"
CALL write_part_mpi(ppos,npartpe,variablefilename,mype)
variablefilename="velocities.bin"
CALL write_part_mpi(pvel,npartpe,variablefilename,mype)
#endif
#ifdef GRAVITY
variablefilename="gravity.bin"
CALL write_var_mpi(phi3d,nx,ny,nz,variablefilename,mype)
#endif

#ifdef NBODY1
do i=1,npartpe,10
   write(UNIT=500+mype,FMT="(4(e13.7,1x))")ppos(1,i),ppos(2,i),ppos(3,i),pvel(1,i)**2+pvel(2,i)**2+pvel(3,i)**2
   write(UNIT=600+mype,FMT="(6(e13.7,1x))")ppos(1,i),ppos(2,i),ppos(3,i),pvel(1,i),pvel(2,i),pvel(3,i)
enddo
#endif
!
! deallocate arrays
!
CALL dealloc_arrays

#ifdef USEMPI
tend = MPI_Wtime()
if(mype == 0)write(*,*)"Total elapsed time: ",tend-tstart
CALL MPI_FINALIZE(ierr)
#endif
!
stop
END PROGRAM hydropad

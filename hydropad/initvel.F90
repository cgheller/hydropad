
SUBROUTINE initvel
!
! initial conditions for the baryonic velocity (linear approximation)
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
!
	IMPLICIT NONE
!
! local variables
#ifdef USEMPI
	INTEGER :: status_array(MPI_STATUS_SIZE,2)
#endif
        INTEGER :: req(2)
        INTEGER :: i,j,k
	INTEGER :: im1,ip1,nxny,to_pe,from_pe
	REAL*8 :: vvv,expon,velocity_factor,z
        REAL*8, DIMENSION (:,:), ALLOCATABLE :: phi_aux
!
	nxny=nx*ny
	allocate(phi_aux(nx,ny))
!
! calculate the correction factor to calculate velocity starting
! from the gravitational force in linear regime
!
	z=redshiftold
	velocity_factor=( (omega_m*(1+z)**3) / &
               (omega_m*(1+z)**3-(omega_m+omega_lambda-1.0)*&
               (1+z)**2+omega_lambda) )**(4.0/7.0)
!
	vvv=0.66666667*velocity_factor
	vvv=vvv*at/dat
	vvv=vvv/omega_m
!
! commented out there is the old version of the initial condition 
! velocity factor
!
!	if(omega0.eq.1.0)then
!	   vvv=sqrt(at)
!	else
!	   expon=omega0**(-0.4)
!	   vvv=(2.0/3.0)*expon/dat
!	endif
!
	   do k=1,nz
	     do j=1,ny
	       do i=1,nx
! vx >
		 ip1=i+1
		 im1=i-1
	         if(i.eq.1)im1=ngrid
	         if(i.eq.ngrid)ip1=1
		 vx3d(i,j,k)=vvv*(phi3d(ip1,j,k)-&
                             phi3d(im1,j,k))*0.5
! vy >
		 ip1=j+1
		 im1=j-1
	         if(j.eq.1)im1=ngrid
	         if(j.eq.ngrid)ip1=1
		 vy3d(i,j,k)=vvv*(phi3d(i,ip1,k)-&
                             phi3d(i,im1,k))*0.5
	       enddo
	     enddo
	   enddo
! vz >
           do k=2,nz-1
             do j=1,ny
               do i=1,nx

		 ip1=k+1
		 im1=k-1
	         if(k.eq.1)im1=ngrid
	         if(k.eq.ngrid)ip1=1
		 vz3d(i,j,k)=vvv*(phi3d(i,j,ip1)-&
                             phi3d(i,j,im1))*0.5
               enddo
             enddo
           enddo
!
#ifdef USEMPI
	   k=1
	   to_pe=mype+1
	   from_pe=mype-1
	   if(mype.eq.npes-1)to_pe=0
	   if(mype.eq.0)from_pe=npes-1
!
	   CALL MPI_Irecv(phi_aux(1,1),nxny,&
                          MPI_DOUBLE_PRECISION,from_pe,10,&
                          MPI_COMM_WORLD,req(1),ierr)
!
	   CALL MPI_Isend(phi3d(1,1,nz),nxny,&
                          MPI_DOUBLE_PRECISION,to_pe,10,MPI_COMM_WORLD,&
                          req(2),ierr)
!
           CALL MPI_WAITALL(2,req,status_array,ierr)
#else

	   phi_aux=phi3d(:,:,nz)   

#endif
!
	   ip1=2
	   do j=1,ny
	     do i=1,nx
	       vz3d(i,j,k)=vvv*(phi3d(i,j,ip1)-&
                           phi_aux(i,j))*0.5
	     enddo
	   enddo
!
#ifdef USEMPI
	   k=nz
	   to_pe=mype-1
	   from_pe=mype+1
	   if(mype.eq.npes-1)from_pe=0
	   if(mype.eq.0)to_pe=npes-1
!
	   CALL MPI_Irecv(phi_aux(1,1),nxny,&
                          MPI_DOUBLE_PRECISION,from_pe,10,&
                          MPI_COMM_WORLD,req(1),ierr)
!
	   CALL MPI_Isend(phi3d(1,1,1),nxny,&
                          MPI_DOUBLE_PRECISION,to_pe,10,MPI_COMM_WORLD,&
                          req(2),ierr)
!
           CALL MPI_WAITALL(2,req,status_array,ierr)
#else

	   phi_aux=phi3d(:,:,1)

#endif
!
	   im1=nz-1
	   do j=1,ny
	     do i=1,nx
	       vz3d(i,j,k)=vvv*(phi_aux(i,j)-&
                           phi3d(i,j,im1))*0.5
	     enddo
	   enddo
!
	deallocate(phi_aux)
! 
END SUBROUTINE initvel

#include "hydrompi.def"

SUBROUTINE dynamics
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
!
! local variables
!
#ifdef USEMPI
        integer :: status_array(MPI_STATUS_SIZE,12)
#endif
	INTEGER :: i,j,k,error
        integer :: req(12)
	integer :: iunit,nunit
	integer :: ipart,i1,i2,i3,ipe3,to_pe,from_pe,countaux
	integer :: n_tot_s,n_tot,nxny
	integer :: index0,stride0
	integer :: index1,stride1
	integer :: i11,i21,i31
	integer :: count_pe
	integer, dimension (:),  allocatable :: count_pe_s
	integer, dimension (:),  allocatable :: count_pe_r
	integer, dimension (:),  allocatable :: point_pe
	integer, dimension (:),  allocatable :: point_pe_r
	integer, dimension (:),  allocatable :: point_pe_aux
	real*8, dimension (:), allocatable :: x1_aux_s
	real*8, dimension (:), allocatable :: x2_aux_s
	real*8, dimension (:), allocatable :: x3_aux_s
	real*8, dimension (:), allocatable :: x1_aux_r
	real*8, dimension (:), allocatable :: x2_aux_r
	real*8, dimension (:), allocatable :: x3_aux_r
	real*8, dimension (:), allocatable :: v1_aux_s
	real*8, dimension (:), allocatable :: v2_aux_s
	real*8, dimension (:), allocatable :: v3_aux_s
	real*8, dimension (:), allocatable :: v1_aux_r
	real*8, dimension (:), allocatable :: v2_aux_r
	real*8, dimension (:), allocatable :: v3_aux_r
	integer, dimension(:), allocatable :: n_pepe
!
	nxny=nx*ny
!
! pre-processing phase
! allocate vectors which depends on the number of processors: npes
!
	allocate (count_pe_s(npes),STAT=error)
	allocate (count_pe_r(npes),STAT=error)
	allocate (point_pe(npes),STAT=error)
	allocate (point_pe_r(npes),STAT=error)
	allocate (point_pe_aux(npes),STAT=error)
        if(error.ne.0)then
          write(*,*)mype,' failed in allocating array in DYNAMICS 1'
          stop
        endif
	count_pe=0
	point_pe=0
	count_pe_s=0
	count_pe_r=0
	point_pe_r=0
!
! find how many particles fall in each processor
!
	n_tot_s=0
	do ipart=1,nparmax
	  i3=int(x3(ipart)+0.5)
          if(i3.eq.0)i3=npes*nz
	  ipe3=(i3-1)/nz
	  if(ipe3.ne.mype)then
	    count_pe_s(ipe3+1)=count_pe_s(ipe3+1)+1
	    n_tot_s=n_tot_s+1
	  endif
	enddo
!
! allocate auxiliary vectors which stores the ids if particles
! belonging to other processors
!
#ifdef DEBUG
	write(*,*)mype,' sends ',n_tot_s,' particles'
#endif

	allocate (x1_aux_s(n_tot_s+1),STAT=error)
	allocate (x2_aux_s(n_tot_s+1),STAT=error)
	allocate (x3_aux_s(n_tot_s+1),STAT=error)
	allocate (v1_aux_s(n_tot_s+1),STAT=error)
	allocate (v2_aux_s(n_tot_s+1),STAT=error)
	allocate (v3_aux_s(n_tot_s+1),STAT=error)
	allocate (n_pepe(n_tot_s),STAT=error)
        if(error.ne.0)then
          write(*,*)mype,' failed in allocating array in DYNAMICS 2'
          stop
        endif
	x1_aux_s=0.0
	x2_aux_s=0.0
	x3_aux_s=0.0
	v1_aux_s=0.0
	v2_aux_s=0.0
	v3_aux_s=0.0
	n_pepe=0
!
! store memory position in auxiliary vectors of the first particle
! which falls in a different pe for each different pe (pointer vectors)
!
	point_pe=0
!	if(mype.ne.0.and.count_pe_s(1).ne.0)point_pe(1)=1
	point_pe(1)=1
!
	countaux=count_pe_s(1)
	do i=2,npes
	      point_pe(i)=point_pe(i-1)+countaux
	      countaux=count_pe_s(i)
	enddo
!
! build auxiliary vectors
!
	point_pe_aux=point_pe
#ifdef USEMPI
	do ipart=1,nparmax
          i3=int(x3(ipart)+0.5)
          if(i3.eq.0)i3=npes*nz
          ipe3=(i3-1)/nz
	  if(ipe3.ne.mype)then
	    index0=point_pe_aux(ipe3+1)
	    x1_aux_s(index0)=x1(ipart)
	    x2_aux_s(index0)=x2(ipart)
	    x3_aux_s(index0)=x3(ipart)
	    v1_aux_s(index0)=v1(ipart)
	    v2_aux_s(index0)=v2(ipart)
	    v3_aux_s(index0)=v3(ipart)
	    n_pepe(index0)=ipart
	    point_pe_aux(ipe3+1)=point_pe_aux(ipe3+1)+1
	  endif
	enddo
!
! communicate pointer vectors
!
	do i=1,npes
	  if(i-1.ne.mype)then
	  to_pe=i-1
	  from_pe=i-1
	  CALL MPI_Send(count_pe_s(i),1,MPI_INTEGER,to_pe,10,&
                        MPI_COMM_WORLD,ierr)
	  CALL MPI_Recv(count_pe_r(i),1,MPI_INTEGER,&
                        from_pe,10,MPI_COMM_WORLD,status,ierr)
	  endif
	enddo
!
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
! calculate the number of particles which falls in each processor
!
	point_pe_r(1)=1

        countaux=count_pe_r(1)
        do i=2,npes
              point_pe_r(i)=point_pe_r(i-1)+countaux
              countaux=count_pe_r(i)
        enddo

	n_tot=0
	do i=1,npes
	  if((i-1).ne.mype)n_tot=n_tot+count_pe_r(i)
	enddo
!
! allocate the destination vectors on each pe
!
#endif
	allocate (x1_aux_r(n_tot+1),STAT=error)
	allocate (x2_aux_r(n_tot+1),STAT=error)
	allocate (x3_aux_r(n_tot+1),STAT=error)
	allocate (v1_aux_r(n_tot+1),STAT=error)
	allocate (v2_aux_r(n_tot+1),STAT=error)
	allocate (v3_aux_r(n_tot+1),STAT=error)
        if(error.ne.0)then
          write(*,*)mype,' failed in allocating array in DYNAMICS 3'
          stop
        endif
	x1_aux_r=0.0
	x2_aux_r=0.0
	x3_aux_r=0.0
	v1_aux_r=0.0
	v2_aux_r=0.0
	v3_aux_r=0.0
!
#ifdef USEMPI
!
! send auxiliary vectors. build destination vectors.
!
	index1=1
	do i=1,npes
	  if((i-1).ne.mype)then 
	    to_pe=i-1
	    from_pe=i-1
	    index0=point_pe(i)
	    index1=point_pe_r(i)
	    stride0=count_pe_s(i)
	    stride1=count_pe_r(i)
!
	    if(stride0.eq.0)then
	      stride0=1
	      index0=n_tot_s+1
	    endif
	    if(stride1.eq.0)then
	      stride1=1
	      index1=n_tot+1
	    endif
!
           CALL MPI_Irecv(x1_aux_r(index1),stride1,MPI_DOUBLE_PRECISION,&
                        from_pe,20,MPI_COMM_WORLD,req(1),ierr)
           CALL MPI_Isend(x1_aux_s(index0),stride0,MPI_DOUBLE_PRECISION,&
                        to_pe,20,MPI_COMM_WORLD,req(2),ierr)

           CALL MPI_Irecv(x2_aux_r(index1),stride1,MPI_DOUBLE_PRECISION,&
                        from_pe,30,MPI_COMM_WORLD,req(3),ierr)
           CALL MPI_Isend(x2_aux_s(index0),stride0,MPI_DOUBLE_PRECISION,&
                        to_pe,30,MPI_COMM_WORLD,req(4),ierr)

           CALL MPI_Irecv(x3_aux_r(index1),stride1,MPI_DOUBLE_PRECISION,&
                        from_pe,40,MPI_COMM_WORLD,req(5),ierr)
           CALL MPI_Isend(x3_aux_s(index0),stride0,MPI_DOUBLE_PRECISION,&
                        to_pe,40,MPI_COMM_WORLD,req(6),ierr)

           CALL MPI_Irecv(v1_aux_r(index1),stride1,MPI_DOUBLE_PRECISION,&
                        from_pe,50,MPI_COMM_WORLD,req(7),ierr)
           CALL MPI_Isend(v1_aux_s(index0),stride0,MPI_DOUBLE_PRECISION,&
                        to_pe,50,MPI_COMM_WORLD,req(8),ierr)

           CALL MPI_Irecv(v2_aux_r(index1),stride1,MPI_DOUBLE_PRECISION,&
                        from_pe,60,MPI_COMM_WORLD,req(9),ierr)
           CALL MPI_Isend(v2_aux_s(index0),stride0,MPI_DOUBLE_PRECISION,&
                        to_pe,60,MPI_COMM_WORLD,req(10),ierr)

           CALL MPI_Irecv(v3_aux_r(index1),stride1,MPI_DOUBLE_PRECISION,&
                        from_pe,70,MPI_COMM_WORLD,req(11),ierr)
           CALL MPI_Isend(v3_aux_s(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,70,MPI_COMM_WORLD,req(12),ierr)
!
            if(npes.gt.1)CALL MPI_WAITALL(12,req,status_array,ierr)
!
	  endif
	enddo
!
#endif
	deallocate (x1_aux_s)
	deallocate (x2_aux_s)
	deallocate (x3_aux_s)
	deallocate (v1_aux_s)
	deallocate (v2_aux_s)
	deallocate (v3_aux_s)
!
	call density(x1_aux_r,x2_aux_r,x3_aux_r,n_tot)
!
	if (nstep.ne.0)then
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! calculate gravitational potential
!
	allocate (gx(nx,ny,0:nz+2),STAT=error)
	allocate (gy(nx,ny,0:nz+2),STAT=error)
	allocate (gz(nx,ny,0:nz+2),STAT=error)
        if(error.ne.0)then
          write(*,*)mype,' failed in allocating array in DYNAMICS 4'
          stop
        endif
!
	call fields
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! calculate gravitational forces
!
	call gravity
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
! move particles
!
	call nbody(n_tot,x1_aux_r,x2_aux_r,x3_aux_r,&
                         v1_aux_r,v2_aux_r,v3_aux_r)
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	deallocate (gx)
	deallocate (gy)
	deallocate (gz)
!
#ifdef USEMPI
	allocate (x1_aux_s(n_tot_s+1),STAT=error)
	allocate (x2_aux_s(n_tot_s+1),STAT=error)
	allocate (x3_aux_s(n_tot_s+1),STAT=error)
	allocate (v1_aux_s(n_tot_s+1),STAT=error)
	allocate (v2_aux_s(n_tot_s+1),STAT=error)
	allocate (v3_aux_s(n_tot_s+1),STAT=error)
        if(error.ne.0)then
          write(*,*)mype,' failed in allocating array in DYNAMICS 5'
          stop
        endif
!
	index0=1
        do i=1,npes
          if((i-1).ne.mype)then
            to_pe=i-1
            from_pe=i-1
            stride0=count_pe_r(i)
            stride1=count_pe_s(i)
	    index1=point_pe(i)
	    index0=point_pe_r(i)
!
	    if(stride0.eq.0)then
	      stride0=1
	      index0=n_tot+1
	    endif
	    if(stride1.eq.0)then
	      stride1=1
	      index1=n_tot_s+1
	    endif
!
           CALL MPI_Irecv(x1_aux_s(index1),stride1,MPI_DOUBLE_PRECISION,&
                         from_pe,20,MPI_COMM_WORLD,req(1),ierr)
           CALL MPI_Isend(x1_aux_r(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,20,MPI_COMM_WORLD,req(2),ierr)

           CALL MPI_Irecv(x2_aux_s(index1),stride1,MPI_DOUBLE_PRECISION,&
                         from_pe,30,MPI_COMM_WORLD,req(3),ierr)
           CALL MPI_Isend(x2_aux_r(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,30,MPI_COMM_WORLD,req(4),ierr)

           CALL MPI_Irecv(x3_aux_s(index1),stride1,MPI_DOUBLE_PRECISION,&
                         from_pe,40,MPI_COMM_WORLD,req(5),ierr)
           CALL MPI_Isend(x3_aux_r(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,40,MPI_COMM_WORLD,req(6),ierr)

           CALL MPI_Irecv(v1_aux_s(index1),stride1,MPI_DOUBLE_PRECISION,&
                         from_pe,50,MPI_COMM_WORLD,req(7),ierr)
           CALL MPI_Isend(v1_aux_r(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,50,MPI_COMM_WORLD,req(8),ierr)

           CALL MPI_Irecv(v2_aux_s(index1),stride1,MPI_DOUBLE_PRECISION,&
                         from_pe,60,MPI_COMM_WORLD,req(9),ierr)
           CALL MPI_Isend(v2_aux_r(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,60,MPI_COMM_WORLD,req(10),ierr)

           CALL MPI_Irecv(v3_aux_s(index1),stride1,MPI_DOUBLE_PRECISION,&
                         from_pe,70,MPI_COMM_WORLD,req(11),ierr)
           CALL MPI_Isend(v3_aux_r(index0),stride0,MPI_DOUBLE_PRECISION,&
                         to_pe,70,MPI_COMM_WORLD,req(12),ierr)
!
            if(npes.gt.1)CALL MPI_WAITALL(12,req,status_array,ierr)
          endif
        enddo
!
	do i=1,n_tot_s
	  x1(n_pepe(i))=x1_aux_s(i)
	  x2(n_pepe(i))=x2_aux_s(i)
	  x3(n_pepe(i))=x3_aux_s(i)
	  v1(n_pepe(i))=v1_aux_s(i)
	  v2(n_pepe(i))=v2_aux_s(i)
	  v3(n_pepe(i))=v3_aux_s(i)
	enddo
!	
#endif
	deallocate (x1_aux_s)
	deallocate (x2_aux_s)
	deallocate (x3_aux_s)
	deallocate (v1_aux_s)
	deallocate (v2_aux_s)
	deallocate (v3_aux_s)
!
	endif
!
#ifdef USEMPI
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
	deallocate (count_pe_s)
	deallocate (count_pe_r)
	deallocate (point_pe)
	deallocate (point_pe_r)
	deallocate (point_pe_aux)
	deallocate (n_pepe)
!	deallocate (x1_aux_s)
!	deallocate (x2_aux_s)
!	deallocate (x3_aux_s)
!	deallocate (v1_aux_s)
!	deallocate (v2_aux_s)
!	deallocate (v3_aux_s)
	deallocate (x1_aux_r)
	deallocate (x2_aux_r)
	deallocate (x3_aux_r)
	deallocate (v1_aux_r)
	deallocate (v2_aux_r)
	deallocate (v3_aux_r)
!
END SUBROUTINE dynamics

#include "hydrompi.def"
!
SUBROUTINE initialdm
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
	INTEGER :: status_array(MPI_STATUS_SIZE,12)
#endif
        INTEGER :: i,j,k
	INTEGER :: irec	
	INTEGER :: nxy
	INTEGER :: to_pe,from_pe
	INTEGER :: req(12)
	REAL*8 :: vxmax,vymax,vzmax,vxmin,vymin,vzmin,vxav,vyav,vzav
	REAL*8 :: xxmax,xymax,xzmax,xxmin,xymin,xzmin
	REAL*8 :: xmpc,vfact
	REAL*8 :: astart
        REAL*8 :: velfactenzo,posfact,convfact
	REAL(KIND=4) ::  x1r4,x2r4,x3r4,v1r4,v2r4,v3r4
        REAL*8, DIMENSION (:), allocatable :: x1_aux_s
        REAL*8, DIMENSION (:), allocatable :: x2_aux_s
        REAL*8, DIMENSION (:), allocatable :: x3_aux_s
        REAL*8, DIMENSION (:), allocatable :: v1_aux_s
        REAL*8, DIMENSION (:), allocatable :: v2_aux_s
        REAL*8, DIMENSION (:), allocatable :: v3_aux_s
!
! read initial expansion factor
!
        enzoactive = 0.0
	nxy=nx*ny
	open(30,file='delta.dat',access='direct',recl=4)
	read(30,rec=1)x1r4
	astart=dble(x1r4)
	read(30,rec=2)x1r4
	velfactenzo=dble(x1r4)
        read(30,rec=3)x1r4
        enzoactive = dble(x1r4)
	xmpc=1.0

        if(enzoactive .eq. -1000.0)then
	   convfact = velfactenzo / 1e5
           posfact = dble(ngrid)
        else
	   convfact = t0h0
           posfact=1.0
        endif

	close(30)
	at=astart
!
! read initial positions and displacements
!
	vxav=0.0
	vyav=0.0
	vzav=0.0
	do i=1,npes
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	if(mype.eq.i-1)then
	  open(30,file='p3m.dat',access='direct',recl=24)
	  do j=1,nparmax
	    irec=j+mype*nparmax
	    read(30,rec=irec)x1r4,x2r4,x3r4,v1r4,v2r4,v3r4
	    x1(j)=dble(x1r4)*posfact
	    x2(j)=dble(x2r4)*posfact
	    x3(j)=dble(x3r4)*posfact
	    v1(j)=dble(v1r4)*convfact
	    v2(j)=dble(v2r4)*convfact
	    v3(j)=dble(v3r4)*convfact
	    vxav=vxav+abs(v1(j))
	    vyav=vyav+abs(v2(j))
	    vzav=vzav+abs(v3(j))
	  enddo
	  close(30)
	endif
	enddo
	vxav=vxav/float(nparmax)
	vyav=vyav/float(nparmax)
	vzav=vzav/float(nparmax)
!
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
        do j=1,nparmax
	  x1(j)=x1(j)*xmpc
	  x2(j)=x2(j)*xmpc
	  x3(j)=x3(j)*xmpc
          if (x1(j).lt.0.0)  x1(j)=x1(j)+xmax
          if (x1(j).ge.xmax) x1(j)=x1(j)-xmax
          if (x2(j).lt.0.0)  x2(j)=x2(j)+xmax
          if (x2(j).ge.xmax) x2(j)=x2(j)-xmax
          if (x3(j).lt.0.0)  x3(j)=x3(j)+xmax
          if (x3(j).ge.xmax) x3(j)=x3(j)-xmax
	enddo


109     format(i2,6(1x,e11.5))


        xxmax=maxval(x1)
        xymax=maxval(x2)
        xzmax=maxval(x3)
        xxmin=minval(x1)
        xymin=minval(x2)
        xzmin=minval(x3)
        vxmax=maxval(v1)
        vymax=maxval(v2)
        vzmax=maxval(v3)
        vxmin=minval(v1)
        vymin=minval(v2)
        vzmin=minval(v3)
	if(mype.eq.0)then
	  write(*,*)
	  write(*,*)'initial data statistics'
	  write(*,*)'maxima, minima of velocity components'
	  write(*,*)'average values of abs velocity components'
	endif
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
        write(*,*)mype,xxmax,xymax,xzmax,xxmin,xymin,xzmin
        write(*,109)mype,vxmax,vymax,vzmax,vxmin,vymin,vzmin
        write(*,109)mype,vxav,vyav,vzav
#ifdef USEMPI
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
	if(npes.gt.1)then
        allocate (x1_aux_s(nxy))
        allocate (x2_aux_s(nxy))
        allocate (x3_aux_s(nxy))
        allocate (v1_aux_s(nxy))
        allocate (v2_aux_s(nxy))
        allocate (v3_aux_s(nxy))
	do j=1,nxy
	  x1_aux_s(j)=x1(j)
	  x2_aux_s(j)=x2(j)
	  x3_aux_s(j)=x3(j)
	  v1_aux_s(j)=v1(j)
	  v2_aux_s(j)=v2(j)
	  v3_aux_s(j)=v3(j)
	enddo
!
	to_pe=mype-1
	if(mype.eq.0)to_pe=npes-1
	from_pe=mype+1
	if(mype.eq.npes-1)from_pe=0
!
        CALL MPI_Irecv(x1(1),nxy,MPI_DOUBLE_PRECISION,&
                      from_pe,10,MPI_COMM_WORLD,req(1),ierr)
        CALL MPI_Isend(x1_aux_s(1),nxy,MPI_DOUBLE_PRECISION,&
                      to_pe,10,MPI_COMM_WORLD,req(2),ierr)

        CALL MPI_Irecv(x2(1),nxy,MPI_DOUBLE_PRECISION,&
                      from_pe,20,MPI_COMM_WORLD,req(3),ierr)
        CALL MPI_Isend(x2_aux_s(1),nxy,MPI_DOUBLE_PRECISION,&
                      to_pe,20,MPI_COMM_WORLD,req(4),ierr)

        CALL MPI_Irecv(x3(1),nxy,MPI_DOUBLE_PRECISION,&
                      from_pe,30,MPI_COMM_WORLD,req(5),ierr)
        CALL MPI_Isend(x3_aux_s(1),nxy,MPI_DOUBLE_PRECISION,&
                      to_pe,30,MPI_COMM_WORLD,req(6),ierr)

        CALL MPI_Irecv(v1(1),nxy,MPI_DOUBLE_PRECISION,&
                      from_pe,40,MPI_COMM_WORLD,req(7),ierr)
        CALL MPI_Isend(v1_aux_s(1),nxy,MPI_DOUBLE_PRECISION,&
                      to_pe,40,MPI_COMM_WORLD,req(8),ierr)

        CALL MPI_Irecv(v2(1),nxy,MPI_DOUBLE_PRECISION,&
                      from_pe,50,MPI_COMM_WORLD,req(9),ierr)
        CALL MPI_Isend(v2_aux_s(1),nxy,MPI_DOUBLE_PRECISION,&
                      to_pe,50,MPI_COMM_WORLD,req(10),ierr)

        CALL MPI_Irecv(v3(1),nxy,MPI_DOUBLE_PRECISION,&
                      from_pe,60,MPI_COMM_WORLD,req(11),ierr)
        CALL MPI_Isend(v3_aux_s(1),nxy,MPI_DOUBLE_PRECISION,&
                      to_pe,60,MPI_COMM_WORLD,req(12),ierr)
!
	CALL MPI_WAITALL(12,req,status_array,ierr)
!
        deallocate (x1_aux_s)
        deallocate (x2_aux_s)
        deallocate (x3_aux_s)
        deallocate (v1_aux_s)
        deallocate (v2_aux_s)
        deallocate (v3_aux_s)
!	  
	endif
#endif
!
END SUBROUTINE initialdm

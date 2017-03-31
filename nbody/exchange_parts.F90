SUBROUTINE exchange_parts(np,ppos,pvel,nrecv,request,nreq)
!
      USE dimension!, ONLY: ngridxpe, ngridype, ngridzpe, nbound
      USE scalar
      USE mpi_inc
      USE nbody_mod
!
! local variables
!
IMPLICIT NONE
!
INTEGER :: np
INTEGER :: i1,i2,i3,iaux
REAL(KIND=8) :: xi1,xi2,xi3
INTEGER :: i,j,k,nreq
INTEGER :: ip1,im1,jp1,j1,jm1,kp1,km1,i3aux,i1aux,i2aux,k1
REAL(KIND=8), INTENT(IN), DIMENSION(3,np) :: ppos
REAL(KIND=8), INTENT(IN), DIMENSION(3,np) :: pvel
REAL(KIND=8) :: x1g,x2g,x3g
INTEGER, DIMENSION(ndims) :: partpe
INTEGER, DIMENSION(ndims) :: boxlaux
INTEGER, DIMENSION(nreq) :: request
INTEGER :: to_pe, from_pe
INTEGER :: nsend,nrecv
INTEGER :: partpeget,reqcount
!
ALLOCATE(nsendpe(npes))
ALLOCATE(nrecvpe(npes))
ALLOCATE(indexsend(npes))
ALLOCATE(psend(npes))
ALLOCATE(precv(npes))
nsendpe=0
nrecvpe=0
indexsend=1
psend=0
precv=0
!
do j=1,np
!
! calculate the cell in which the particle lies in local box coordinates
!
  partpe=0
  partpeget = 0
#ifdef USEMPI
  partpe(1) = int(ppos(1,j)*npesx)
  partpe(2) = int(ppos(2,j)*npesy)
  partpe(3) = int(ppos(3,j)*npesz)
  CALL MPI_CART_RANK(COMM_CART, partpe, partpeget, ierr)
#endif
  if(partpeget .ne. mype)nsendpe(partpeget) = nsendpe(partpeget)+1
enddo
!
! distribute information on particles to exchange
!
CALL MPI_ALLTOALL(nsendpe,1,MPI_INTEGER,nrecvpe,1,MPI_INTEGER,COMM_CART,ierr)
!
! allocate arrays for exchange particles data
!
nsend = 0
nrecv = 0
do i=1,npes
  nsend = nsend+nsendpe(i)
  nrecv = nrecv+nrecvpe(i)
enddo 
ALLOCATE(xsend(3,nsend))
ALLOCATE(vsend(3,nsend))
ALLOCATE(xrecv(3,nrecv))
ALLOCATE(vrecv(3,nrecv))
do i=1,npes-1
  psend(i+1) = psend(i)+nsendpe(i)
  precv(i+1) = precv(i)+nrecvpe(i)
enddo 
! 
! build the arrays to communicate
!
do j=1,np
!
! calculate the cell in which the particle lies in local box coordinates
!
  partpe=0
  partpeget = 0
#ifdef USEMPI
  partpe(1) = int(ppos(1,j)*npesx)
  partpe(2) = int(ppos(2,j)*npesy)
  partpe(3) = int(ppos(3,j)*npesz)
  CALL MPI_CART_RANK(COMM_CART, partpe, partpeget, ierr)
#endif
!
  if(partpeget == mype)CYCLE
!
! fill arrays
!
  iaux = psend(partpeget)+indexsend(partpeget)
  indexsend(partpeget) = indexsend(partpeget)+1
  xsend(1,iaux)=ppos(1,j)
  xsend(2,iaux)=ppos(2,j)
  xsend(3,iaux)=ppos(3,j)
  vsend(1,iaux)=pvel(1,j)
  vsend(2,iaux)=pvel(2,j)
  vsend(3,iaux)=pvel(3,j)
!
enddo
!
! Communication
!
reqcount=1
do i=1,npes
  to_pe = i-1
  from_pe = i-1
  CALL MPI_ISEND(xsend(1,psend(i)+1),3*nsendpe(i),MPI_DOUBLE_PRECISION,&
                 to_pe,10,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
  CALL MPI_IRECV(xrecv(1,precv(i)+1),3*nrecvpe(i),MPI_DOUBLE_PRECISION,&
                 from_pe,10,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
  CALL MPI_ISEND(vsend(1,psend(i)+1),3*nsendpe(i),MPI_DOUBLE_PRECISION,&
                 to_pe,20,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
  CALL MPI_IRECV(vrecv(1,precv(i)+1),3*nrecvpe(i),MPI_DOUBLE_PRECISION,&
                 from_pe,20,COMM_CART,request(reqcount),ierr)
enddo
!
END SUBROUTINE exchange_parts

SUBROUTINE exchangeback_parts(np,ppos,pvel)
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
INTEGER :: i,j,k,nreq
REAL(KIND=8), INTENT(INOUT), DIMENSION(3,np) :: ppos
REAL(KIND=8), INTENT(INOUT), DIMENSION(3,np) :: pvel
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: statcomm
INTEGER, DIMENSION(ndims) :: partpe
INTEGER, DIMENSION(ndims) :: boxlaux
INTEGER, ALLOCATABLE, DIMENSION(:) :: request
INTEGER :: to_pe, from_pe
INTEGER :: partpeget,reqcount,partpeg
!
nreq = 4*npes
ALLOCATE(request(nreq))
ALLOCATE(statcomm(MPI_STATUS_SIZE,nreq))
!
! Communication: send back updated data
!
reqcount=1
do i=1,npes
  to_pe = i-1
  from_pe = i-1
  CALL MPI_ISEND(xrecv(1,precv(i)+1),3*nrecvpe(i),MPI_DOUBLE_PRECISION,&
                 to_pe,10,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
  CALL MPI_IRECV(xsend(1,psend(i)+1),3*nsendpe(i),MPI_DOUBLE_PRECISION,&
                 from_pe,10,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
  CALL MPI_ISEND(vrecv(1,precv(i)+1),3*nrecvpe(i),MPI_DOUBLE_PRECISION,&
                 to_pe,20,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
  CALL MPI_IRECV(vsend(1,psend(i)+1),3*nsendpe(i),MPI_DOUBLE_PRECISION,&
                 from_pe,20,COMM_CART,request(reqcount),ierr)
  reqcount=reqcount+1
enddo
!
! complete Nbody communication
!
CALL MPI_Waitall(nreq,request,statcomm,ierr)
!
! restore data in original locations
!
iaux=1
do j=1,nsend
!
! fill arrays
!
  ppos(1,listofparticles(j))=xsend(1,j)
  ppos(2,listofparticles(j))=xsend(2,j)
  ppos(3,listofparticles(j))=xsend(3,j)
  pvel(1,listofparticles(j))=vsend(1,j)
  pvel(2,listofparticles(j))=vsend(2,j)
  pvel(3,listofparticles(j))=vsend(3,j)
!
enddo
!
DEALLOCATE(request)
DEALLOCATE(statcomm)
!
END SUBROUTINE exchangeback_parts

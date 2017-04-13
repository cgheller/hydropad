SUBROUTINE drift(np,ppos,pvel)
!
      USE dimension
      USE scalar
      USE mpi_inc
!
! local variables
!
IMPLICIT NONE
!
INTEGER :: i,j,k,np
REAL(KIND=8), INTENT(INOUT), DIMENSION(3,np) :: ppos
REAL(KIND=8), INTENT(INOUT), DIMENSION(3,np) :: pvel
!
! move particles using a 3 steps "KDK" method: 
!                 1)  v^{n+1/2}   = v^{n} + h/2 * F^{n}
!                 2)  x^{n+1}     = x^{n} + h v^{n+1/2}
!                 3)  v^{n+1}     = v^{n+1/2} + h/2 * F^{n+1}
! STEP 2: particles drift
!
do j=1,np
!
          if(ppos(1,j) == -1)CYCLE
          ppos(1,j)=ppos(1,j)+rdtath*pvel(1,j)
          ppos(2,j)=ppos(2,j)+rdtath*pvel(2,j)
          ppos(3,j)=ppos(3,j)+rdtath*pvel(3,j)
          if (ppos(1,j).lt.0.0) ppos(1,j)=ppos(1,j)+boxsize(1)
          if (ppos(1,j).ge.boxsize(1)) ppos(1,j)=ppos(1,j)-boxsize(1)
          if (ppos(2,j).lt.0.0) ppos(2,j)=ppos(2,j)+boxsize(2)
          if (ppos(2,j).ge.boxsize(2)) ppos(2,j)=ppos(2,j)-boxsize(2)
          if (ppos(3,j).lt.0.0) ppos(3,j)=ppos(3,j)+boxsize(3)
          if (ppos(3,j).ge.boxsize(3)) ppos(3,j)=ppos(3,j)-boxsize(3)
!
enddo
!
END SUBROUTINE drift

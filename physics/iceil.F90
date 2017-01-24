INTEGER FUNCTION iceil(i,j)

      IMPLICIT NONE

      INTEGER :: i,j
      REAL*8  :: a

      a = dble(i)/dble(j)

      iceil = CEILING(a)

END FUNCTION iceil


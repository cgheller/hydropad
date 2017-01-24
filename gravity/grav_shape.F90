SUBROUTINE grav_shape
!
!  shape calculates the shape function used in the Green's function potential
!  calculation.  This function smooths the force law to give the optimal
!  force for a given particle shape using exact Fourier space gradient and
!  CIC interpolation.  See Hockney and Eastwood, section 8-3.
!  2*lmax+1 gives the number of aliases to include in the shape function.
!  eta is the grid force softening length, in units of the grid spacing.
!
      USE dimension
      USE matrix
      USE vector
      USE scalar
!
      USE mpi_inc
!
! local variables
!
      INTEGER :: i,j,k
      INTEGER :: iunit
      INTEGER :: ishape,nshap,lmax,k3,k3aux,k2,k1,iindex,l3,l2,l1
      REAL*8 :: c23,eta,etahalf,arg1,arg0,arg,sinarg,a,a3,a2,a1
      REAL*8 :: sum1,sum2,sum3,b3,b33,a3b3,b2,b22,b23,a2b2,b1,b11
      REAL*8 :: b123,a1b1,qq,shat,q
      REAL*8 :: utemp(2048)
      CHARACTER*50 :: shapename
      PARAMETER (c23=2.0e0/3.0e0)
!
      nshap=ngrid
      lmax=2
      eta=3.5
!
      etahalf=0.5e0*eta
      arg1=twopi/ngrid
      arg0=pi/ngrid
      utemp(1)=1.0e0
      do 10 k=1,n11
!
!  Calculate normalized 1-D CIC assignment function.
!
      arg=k*arg0
      if (k.gt.n12) arg=arg-pi
      sinarg=sin(arg)
      a=sinarg/arg
      utemp(k+1)=a*a/(1.0e0-c23*sinarg*sinarg)
10    continue
!
!  Calculate "influence" function for particle with linear density profile
!  (Hockney and Eastwood shape S2(r)).
!
      do 70 k3=0,n11pe
      k3aux = k3 + mype*n1pe
      a3=k3aux*arg1
      if (k3aux.ge.n12) a3=a3-twopi
        do 60 k2=0,n11
        a2=k2*arg1
        if (k2.ge.n12) a2=a2-twopi
          do 50 k1=0,n11
          a1=k1*arg1
          if (k1.ge.n12) a1=a1-twopi
          iindex=1+k1+k2*n1+k3*n1*n1
          if (iindex.eq.1.and.mype.eq.0) then
            grav_shap(k1,k2,k3)=1.0
            go to 50
          end if
          a=utemp(k1+1)*utemp(k2+1)*utemp(k3aux+1)
          grav_shap(k1,k2,k3)=a*a
          sum1=0.0e0
          sum2=0.0e0
          sum3=0.0e0
            do 40 l3=-lmax,lmax
            b3=a3+l3*twopi
            b33=b3*b3
            a3b3=1.0e0
            if (l3.ne.0) a3b3=a3/b3
              do 30 l2=-lmax,lmax
              b2=a2+l2*twopi
              b22=b2*b2
              b23=b22+b33
              a2b2=1.0e0
              if (l2.ne.0) a2b2=a2/b2
                do 20 l1=-lmax,lmax
                b1=a1+l1*twopi
                b11=b1*b1
                b123=b11+b23
                a1b1=1.0e0
                if (l1.ne.0) a1b1=a1/b1
                q=etahalf*sqrt(b123)
                qq=q*q
                shat=12.0e0/(qq*qq)*(2.0e0-2.0e0*cos(q)-q*sin(q))
                a=a1b1*a2b2*a3b3
                a=a*a*shat
                a=a*a/b123
                sum1=sum1+b1*a
                sum2=sum2+b2*a
                sum3=sum3+b3*a
20              continue
30            continue
40          continue
          grav_shap(k1,k2,k3)=grav_shap(k1,k2,k3)*&
                              (a1*sum1+a2*sum2+a3*sum3)
50        continue
60      continue
70    continue
!
END SUBROUTINE grav_shape

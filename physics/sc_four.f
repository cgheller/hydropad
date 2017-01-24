      subroutine fft3rinv(a,n1)
c  fft3rinv performs a 3-D inverse FFT from the spatial domain to the spectral
c  domain.   Array a is real of length n1*n1*(n1+2).
c  On input, the first n1*n1*n1 elements are filled with values of a(i,j,k)
c  in the spatial domain.  On output, the first n1*n1*n1 elements are filled
c  with values a(i,j,k), i=1,...,n1/2 of the complex Fourier transform.
c  The last 2*n1*n1 values are filled with complex transform values for
c  i=n/2+1.  fft3rinv is based on FOURN from Press et al., Numerical Recipes.
c  N.B. The transform from spatial to spectral domain is defined with a
c  minus sign in the complex exponential.  To change this, change the sign
c  of twopi.  After transforming and then performing the inverse transform,
c  a(i,j,k) must be divided by n1*n1*n1.
      implicit real*4 (a-h,o-z)
      parameter (twopi=-6.283185307179586d0)
      dimension a(*)
      real*8 theta,wpr,wpi,wr,wi,wtemp
c
      ntot=n1*n1*n1/2
      nprev=1
c  Main loop over the dimensions.
      do 18 idim=1,3
      n=n1
      if (idim.eq.1) n=n1/2
      nrem=ntot/(n*nprev)
      ip1=2*nprev
      ip2=ip1*n
      ip3=ip2*nrem
      i2rev=1
c  This is the bit reversal section of the routine.
        do 14 i2=1,ip2,ip1
        if (i2.lt.i2rev) then
          do 13 i1=i2,i2+ip1-2,2
            do 12 i3=i1,ip3,ip2
            i3rev=i2rev+i3-i2
            tempr=a(i3)
            tempi=a(i3+1)
            a(i3)=a(i3rev)
            a(i3+1)=a(i3rev+1)
            a(i3rev)=tempr
            a(i3rev+1)=tempi
12          continue
13        continue
        end if
        ibit=ip2/2
1        if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
          i2rev=i2rev-ibit
          ibit=ibit/2
          go to 1
        end if
        i2rev=i2rev+ibit
14      continue
c  Here begins the Danielson-Lanczos section of the routine.
      ifp1=ip1
2      if (ifp1.lt.ip2) then
        ifp2=2*ifp1
c  Initialize for the trigonometric recurrence.
        theta=twopi/(ifp2/ip1)
        wpr=sin(0.5d0*theta)
        wpr=-2.0d0*wpr*wpr
        wpi=sin(theta)
        wr=1.0d0
        wi=0.0d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
c  Danielson-Lanczos formula.
              k1=i2
              k2=k1+ifp1
              tempr=wr*a(k2)-wi*a(k2+1)
              tempi=wr*a(k2+1)+wi*a(k2)
              a(k2)=a(k1)-tempr
              a(k2+1)=a(k1+1)-tempi
              a(k1)=a(k1)+tempr
              a(k1+1)=a(k1+1)+tempi
15            continue
16          continue
c  Trigonometric recurrence.
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
17        continue
        ifp1=ifp2
        go to 2
      end if
      nprev=n*nprev
18    continue
c  Use symmetries to obtain the true transform.  This part is specially
c  written for 3-D.
c  Initialize for the trigonometric recurrence.
      theta=twopi/n1
      wpr=sin(0.5d0*theta)
      wpr=-2.0d0*wpr*wpr
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
c  Save j1=1,2 for later.
      np2=n1+2
      n21=np2/2
      n12=n1*n1
      n13=n12*n1
      do 30 j1=3,n21,2
      k1=np2-j1
      i1=j1
      i0=k1+n12+n13
        do 28 j3=1,n1
          do 26 j2=1,n1
          i2=i0
          if (j2.eq.1) i2=i2-n12
          if (j3.eq.1) i2=i2-n13
          if ((j1.eq.k1).and.(i1.gt.i2)) go to 25
          ffer=0.5*(a(i1)+a(i2))
          ffei=0.5*(a(i1+1)-a(i2+1))
          ffor=0.5*(a(i1+1)+a(i2+1))
          ffoi=-0.5*(a(i1)-a(i2))
          tempr=wr*ffor-wi*ffoi
          tempi=wi*ffor+wr*ffoi
          a(i1)=ffer+tempr
          a(i1+1)=ffei+tempi
          a(i2)=ffer-tempr
          a(i2+1)=-ffei+tempi
25          i1=i1+n1
          i0=i0-n1
26        continue
28      continue
c  Trigonometric recurrence.
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
30    continue
c  j1=1,2.
      ij1=0
      i0=n1+n12
      do 34 j3=1,n1
        do 32 j2=1,n1
        ij2=i0
        if (j2.eq.1) ij2=ij2-n1
        if (j3.eq.1) ij2=ij2-n12
        if (ij1.gt.ij2) go to 31
        i1=1+n1*ij1
        i2=1+n1*ij2
        i3=n13+1+2*ij1
        i4=n13+1+2*ij2
        ffer=0.5*(a(i1)+a(i2))
        ffei=0.5*(a(i1+1)-a(i2+1))
        ffor=0.5*(a(i1+1)+a(i2+1))
        ffoi=-0.5*(a(i1)-a(i2))
        a(i1)=ffer+ffor
        a(i1+1)=ffei+ffoi
        a(i2)=a(i1)
        a(i2+1)=-a(i1+1)
        a(i4)=ffer-ffor
        a(i4+1)=-ffei+ffoi
        a(i3)=a(i4)
        a(i3+1)=-a(i4+1)
31        ij1=ij1+1
        i0=i0-1
32        continue
34      continue
      return
      end
C
C*****************************************************************
C*****************************************************************
C
      subroutine fft3r(a,n1)
c  fft3r performs a 3-D FFT from the spectral domain to the spatial domain.
c  Array a is real of length n1*n1*(n1+2).  On input, the first n1*n1*n1
c  elements are filled with values a(i,j,k), i=1,...,n1/2 of the complex
c  Fourier transform.  The last 2*n1*n1 values are filled with complex
c  transform values for i=n1/2+1.  On output, the first n1*n1*n1 elements are
c  filled with values of a(i,j,k) in the spatial domain.
c  fft3r is based on FOURN from Press et al., Numerical Recipes.
c  N.B. The transform from spectral to spatial domain is defined with a
c  plus sign in the complex exponential.  To change this, change the sign
c  of twopi.  After transforming and then performing the inverse transform,
c  a(i,j,k) must be divided by n1*n1*n1.
      implicit real*4 (a-h,o-z)
      parameter (twopi=6.283185307179586d0)
      dimension a(*)
      real*8 theta,wpr,wpi,wr,wi,wtemp
c
c  Use symmetries to obtain the true transform.  This part is specially
c  written for 3-D.
c  Initialize for the trigonometric recurrence.
      theta=twopi/n1
      wpr=sin(0.5d0*theta)
      wpr=-2.0d0*wpr*wpr
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
c  Save j1=1,2 for later.
      np2=n1+2
      n21=np2/2
      n12=n1*n1
      n13=n12*n1
      do 30 j1=3,n21,2
      k1=np2-j1
      i1=j1
      i0=k1+n12+n13
        do 28 j3=1,n1
          do 26 j2=1,n1
          i2=i0
          if (j2.eq.1) i2=i2-n12
          if (j3.eq.1) i2=i2-n13
          if ((j1.eq.k1).and.(i1.gt.i2)) go to 25
c  Multiply by 2 since this is a real FFT summing over only half
c  of the harmonics.
          ffer=a(i1)+a(i2)
          ffei=a(i1+1)-a(i2+1)
          ffor=-(a(i1+1)+a(i2+1))
          ffoi=a(i1)-a(i2)
          tempr=wr*ffor-wi*ffoi
          tempi=wi*ffor+wr*ffoi
          a(i1)=ffer+tempr
          a(i1+1)=ffei+tempi
          a(i2)=ffer-tempr
          a(i2+1)=-ffei+tempi
25          i1=i1+n1
          i0=i0-n1
26        continue
28      continue
c  Trigonometric recurrence.
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
30    continue
c  j1=1,2.
      ij1=0
      i0=n1+n12
      do 34 j3=1,n1
        do 32 j2=1,n1
        ij2=i0
        if (j2.eq.1) ij2=ij2-n1
        if (j3.eq.1) ij2=ij2-n12
        if (ij1.gt.ij2) go to 31
        i1=1+n1*ij1
        i2=1+n1*ij2
        i4=n13+1+2*ij2
        ffer=a(i1)+a(i4)
        ffei=a(i1+1)-a(i4+1)
        ffor=-(a(i1+1)+a(i4+1))
        ffoi=a(i1)-a(i4)
        a(i1)=ffer+ffor
        a(i1+1)=ffei+ffoi
        a(i2)=ffer-ffor
        a(i2+1)=-ffei+ffoi
31        ij1=ij1+1
        i0=i0-1
32        continue
34      continue
c
      ntot=n1*n1*n1/2
      nprev=1
c  Main loop over the dimensions.
      do 18 idim=1,3
      n=n1
      if (idim.eq.1) n=n1/2
      nrem=ntot/(n*nprev)
      ip1=2*nprev
      ip2=ip1*n
      ip3=ip2*nrem
      i2rev=1
c  This is the bit reversal section of the routine.
        do 14 i2=1,ip2,ip1
        if (i2.lt.i2rev) then
          do 13 i1=i2,i2+ip1-2,2
            do 12 i3=i1,ip3,ip2
            i3rev=i2rev+i3-i2
            tempr=a(i3)
            tempi=a(i3+1)
            a(i3)=a(i3rev)
            a(i3+1)=a(i3rev+1)
            a(i3rev)=tempr
            a(i3rev+1)=tempi
12          continue
13        continue
        end if
        ibit=ip2/2
1        if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
          i2rev=i2rev-ibit
          ibit=ibit/2
          go to 1
        end if
        i2rev=i2rev+ibit
14      continue
c  Here begins the Danielson-Lanczos section of the routine.
      ifp1=ip1
2      if (ifp1.lt.ip2) then
        ifp2=2*ifp1
c  Initialize for the trigonometric recurrence.
        theta=twopi/(ifp2/ip1)
        wpr=sin(0.5d0*theta)
        wpr=-2.0d0*wpr*wpr
        wpi=sin(theta)
        wr=1.0d0
        wi=0.0d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
c  Danielson-Lanczos formula.
              k1=i2
              k2=k1+ifp1
              tempr=wr*a(k2)-wi*a(k2+1)
              tempi=wr*a(k2+1)+wi*a(k2)
              a(k2)=a(k1)-tempr
              a(k2+1)=a(k1+1)-tempi
              a(k1)=a(k1)+tempr
              a(k1+1)=a(k1+1)+tempi
15            continue
16          continue
c  Trigonometric recurrence.
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
17        continue
        ifp1=ifp2
        go to 2
      end if
      nprev=n*nprev
18    continue
      return
      end
C

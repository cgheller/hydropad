SUBROUTINE gravfield(A,iphig)
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE mpi_inc
!
! local variables
!
	IMPLICIT NONE
!
        INTEGER :: i,j,k,error
	INTEGER :: iphig
	INTEGER :: k_pe,kaux,kshift
	INTEGER :: ip1,jp1,kp1
	REAL*8  :: ga,akk,ainvakk
        COMPLEX*16 :: A(nx,ny,nz)
        COMPLEX*16, DIMENSION (:,:,:), ALLOCATABLE :: B
!
        allocate(B(nx,ny,nz),STAT=error)
        if(error.ne.0)then
          write(*,*)mype,' failed in allocating array B in GRAVFIELDS'
          stop
        endif
!
	if(iphig.eq.0)then
	  A=dcmplx(phi3d,0.0)
	  call pc3fft(A,desca,1)
          if (mype.eq.0) then
	     A(1,1,1) = (0.0, 0.0)
	     B(1,1,1) = (0.0, 0.0)
	  endif
!
! calculate gravitational field >
!
          k_pe = mype * n1pe
          ga=gconst/at
          do k = 0, nz-1
	    kp1=k+1
            do j = 0, ny-1
	      jp1=j+1
              do i = 1, nx-1
	        ip1=i+1
                kaux = k + k_pe
                akk = akq(i) + akq(j) + akq(kaux)
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = A(ip1,jp1,kp1) * ga *&
                                 grav_shap(i,j,k) * ainvakk
              enddo
            enddo
          enddo
!
! fa il piano x = n12 = n1/2
!
!        do k = 0, n11pe
!	  kp1=k+1
!          do j = 0, n11
!	    jp1=j+1
!            i = n12
!	    ip1=i+1
!            kaux = k + k_pe
!            akk = akq(i) + akq(j) + akq(kaux)
!            ainvakk = 1. / akk
!            B(ip1,jp1,kp1) =  A(ip1,jp1,kp1) * ga *&
!                             grav_shap(i,j,k) * ainvakk
!          enddo
!        enddo
!
! fa il piano x=0
!
          do k = 0, n11pe
	    kp1=k+1
            do j = 0, n11
	      jp1=j+1
              i = 0
	      ip1=i+1
              kaux = k + k_pe
              akk = akq(i) + akq(j) + akq(kaux)
              if (akk .ne. 0.0) then
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = A(ip1,jp1,kp1) * ga *&
                                 grav_shap(i,j,k) * ainvakk
              endif
            enddo
          enddo
!
! FFT transform from k to real-space
!
	  call pc3fft(B,desca,-1)
!
	  phi3d=real(B)
!	  phi3d=dreal(B)
!
! calculate gravitational forces : x component >
!
	else if(iphig.eq.1)then
          if (mype.eq.0) B(1,1,1) = (0.0, 0.0)
          k_pe = mype * n1pe
          ga=gconst/at
          do k = 0, nz-1
            kp1=k+1
            do j = 0, ny-1
              jp1=j+1
              do i = 1, nx-1
                ip1=i+1
                kaux = k + k_pe
                akk = akq(i) + akq(j) + akq(kaux)
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = ak(i) * dcmplx(0.0,1.0) *& 
                                 A(ip1,jp1,kp1) *&
                                 ga * grav_shap(i,j,k) * ainvakk
              enddo
            enddo
          enddo
!
          do k = 0, n11pe
            kp1=k+1
            do j = 0, n11
              jp1=j+1
              i = 0
              ip1=i+1
              kaux = k + k_pe
              akk = akq(i) + akq(j) + akq(kaux)
              if (akk .ne. 0.0) then
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = ak(i) * dcmplx(0.0,1.0) *& 
                                 A(ip1,jp1,kp1) *&
                                 ga * grav_shap(i,j,k) * ainvakk
              endif
            enddo
          enddo
!
	call pc3fft(B,desca,-1)
!
	do k=1,nz
	do j=1,ny
	do i=1,nx
	  gx(i,j,k)=-real(B(i,j,k))
!	  gx(i,j,k)=-dreal(B(i,j,k))
	enddo
	enddo
	enddo
!
! calculate gravitational forces : y component >
!
	else if(iphig.eq.2)then
          if (mype.eq.0) B(1,1,1) = (0.0, 0.0)
          k_pe = mype * n1pe
          ga=gconst/at
          do k = 0, nz-1
            kp1=k+1
            do j = 0, ny-1
              jp1=j+1
              do i = 1, nx-1
                ip1=i+1
                kaux = k + k_pe
                akk = akq(i) + akq(j) + akq(kaux)
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = ak(j) * dcmplx(0.0,1.0) *& 
                                 A(ip1,jp1,kp1) *&
                                 ga * grav_shap(i,j,k) * ainvakk
              enddo
            enddo
          enddo
!
          do k = 0, n11pe
            kp1=k+1
            do j = 0, n11
              jp1=j+1
              i = 0
              ip1=i+1
              kaux = k + k_pe
              akk = akq(i) + akq(j) + akq(kaux)
              if (akk .ne. 0.0) then
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = ak(j) * dcmplx(0.0,1.0) * &
                                 A(ip1,jp1,kp1) *&
                                 ga * grav_shap(i,j,k) * ainvakk
              endif
            enddo
          enddo
!
	call pc3fft(B,desca,-1)
!
	do k=1,nz
	do j=1,ny
	do i=1,nx
	  gy(i,j,k)=-real(B(i,j,k))
!	  gy(i,j,k)=-dreal(B(i,j,k))
	enddo
	enddo
	enddo
!
! calculate gravitational forces : z component >
!
	else if(iphig.eq.3)then
          if (mype.eq.0) B(1,1,1) = (0.0, 0.0)
          k_pe = mype * n1pe
          ga=gconst/at
          do k = 0, nz-1
            kp1=k+1
            do j = 0, ny-1
              jp1=j+1
              do i = 1, nx-1
                ip1=i+1
                kaux = k + k_pe
                akk = akq(i) + akq(j) + akq(kaux)
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = ak(kaux) * dcmplx(0.0,1.0) *& 
                                 A(ip1,jp1,kp1) *&
                                 ga * grav_shap(i,j,k) * ainvakk
              enddo
            enddo
          enddo
!
          do k = 0, n11pe
            kp1=k+1
            do j = 0, n11
              jp1=j+1
              i = 0
              ip1=i+1
              kaux = k + k_pe
              akk = akq(i) + akq(j) + akq(kaux)
              if (akk .ne. 0.0) then
                ainvakk = 1. / akk
                B(ip1,jp1,kp1) = ak(kaux) * dcmplx(0.0,1.0) *& 
                                 A(ip1,jp1,kp1) *&
                                 ga * grav_shap(i,j,k) * ainvakk
              endif
            enddo
          enddo
!
	call pc3fft(B,desca,-1)
!
	do k=1,nz
	do j=1,ny
	do i=1,nx
	  gz(i,j,k)=-real(B(i,j,k))
!	  gz(i,j,k)=-dreal(B(i,j,k))
	enddo
	enddo
	enddo
!
        endif
!
	deallocate(B)
!
END SUBROUTINE gravfield

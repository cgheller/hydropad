SUBROUTINE gforce
!
! This routine compute the value of the components
! of the gravitational force starting from the potential
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
!
! local variables
!
	IMPLICIT NONE
        INTEGER :: i,j,k
	REAL*8 :: phihalf(ngrid)
	REAL*8 :: drho,dt2dt
!
	dt2dt=dt/(2.0*dtold)
!
! calculate gravitational potential at half time step 
!
 	do i=1,ngrid
	  phihalf(i)=phi0(i)
	enddo
!
! compute gravitational force
!
	   i=1
	    if((rho(i)-rho(ngrid))*&
               (rho(i+1)-rho(i))&
               .gt.0.0.and.(rho(i+1)-&
               rho(ngrid)).ne.0.0)then
	       drho=min(&
               abs(rho(i+1)-rho(ngrid)),&
               2.0*abs(rho(i)-rho(ngrid)),&
               2.0*abs(rho(i+1)-rho(i)))*&
               (rho(i+1)-rho(ngrid))/&
               abs(rho(i+1)-rho(ngrid))
	     else
	      drho=0.0
	     endif
	     g(i)=(0.5*rdx)*&
                  (phi1(i+1)-phi1(ngrid)+&
                  (phi1(i+1)-2.0*phi1(i)+&
                   phi1(ngrid))*drho/(12.0*&
                   rho(i)))
	     ghalf(i)=(0.5*rdx)*&
                  (phihalf(i+1)-phihalf(ngrid)+&
                  (phihalf(i+1)-2.0*phihalf(i)+&
                   phihalf(ngrid))*drho/(12.0*&
                   rho(i)))
	   i=ngrid
	    if((rho(i)-rho(i-1))*&
               (rho(1)-rho(i))&
               .gt.0.0.and.(rho(1)-rho(i-1))&
               .ne.0.0)then
	       drho=min(&
               abs(rho(1)-rho(i-1)),&
               2.0*abs(rho(i)-rho(i-1)),&
               2.0*abs(rho(1)-rho(i)))*&
               (rho(1)-rho(i-1))/&
               abs(rho(1)-rho(i-1))
	     else
	      drho=0.0
	     endif
	     g(i)=(0.5*rdx)*&
                  (phi1(1)-phi1(i-1)+&
                  (phi1(1)-2.0*phi1(i)+&
                   phi1(i-1))*drho/(12.0*&
                   rho(i)))
	     ghalf(i)=(0.5*rdx)*&
                  (phihalf(1)-phihalf(i-1)+&
                  (phihalf(1)-2.0*phihalf(i)+&
                   phihalf(i-1))*drho/(12.0*&
                   rho(i)))
	   do i=2,ngrid-1
	    if((rho(i)-rho(i-1))*&
               (rho(i+1)-rho(i))&
               .gt.0.0.and.(rho(i+1)-&
               rho(i-1)).ne.0.0)then
	       drho=min(&
               abs(rho(i+1)-rho(i-1)),&
               2.0*abs(rho(i)-rho(i-1)),&
               2.0*abs(rho(i+1)-rho(i)))*&
               (rho(i+1)-rho(i-1))/&
               abs(rho(i+1)-rho(i-1))
	     else
	      drho=0.0
	     endif
	     g(i)=(0.5*rdx)*&
                  (phi1(i+1)-phi1(i-1)+&
                  (phi1(i+1)-2.0*phi1(i)+&
                   phi1(i-1))*drho/(12.0*&
                   rho(i)))
	     ghalf(i)=(0.5*rdx)*&
                  (phihalf(i+1)-phihalf(i-1)+&
                  (phihalf(i+1)-2.0*phihalf(i)+&
                   phihalf(i-1))*drho/(12.0*&
                   rho(i)))
	   enddo
!
END SUBROUTINE gforce

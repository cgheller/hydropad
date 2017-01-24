SUBROUTINE nbody(n_tot,x1_aux_r,x2_aux_r,x3_aux_r,v1_aux_r,v2_aux_r,v3_aux_r)
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
      INTEGER :: i,j,k
      INTEGER :: n_tot
      INTEGER :: ip1,im1,jp1,j1,jm1,kp1,km1,i3aux,i1aux,i2aux
      INTEGER :: iidim1,mype1
      INTEGER :: ip,jp,kp,ind1p,ind2p,ind0p1,ind0p2,ind0m1
      INTEGER :: i11,i21,i31,i1,i2,i3,ipe3
      INTEGER :: ipp1,ipm1
      INTEGER :: ind,ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8
      INTEGER :: k1,k2,k3,k23,k_pe
!
      REAL*8  :: dk,akmax,sigma,dddmax,dtat
      REAL*8  :: dd1,dd2,dd3,de1,de2,de3
      REAL*8  :: ak3,ak33,ak2,ak23,ak1,akk,ffx,ffy,ffz,rhommm
      REAL*8  :: invngrid3
      REAL*8  :: xhalf(3),vhalf(3),force(8)
      REAL*8  :: dephi,phih
      REAL*8  :: rhotot,dm_mean_tot,bm_mean_tot,dm_mean,bm_mean
      REAL*8  :: de10,de20,de30,dd10,dd20,dd30
      REAL*8  :: d1,d2,d3,d4,d5,d6,d7,d8
!
      REAL*8, INTENT(INOUT) :: x1_aux_r(n_tot+1)
      REAL*8, INTENT(INOUT) :: x2_aux_r(n_tot+1)
      REAL*8, INTENT(INOUT) :: x3_aux_r(n_tot+1)
      REAL*8, INTENT(INOUT) :: v1_aux_r(n_tot+1)
      REAL*8, INTENT(INOUT) :: v2_aux_r(n_tot+1)
      REAL*8, INTENT(INOUT) :: v3_aux_r(n_tot+1)
!
      dk=twopi/xmax
      akmax=ngrid*dk
      mype1=mype+1
      invngrid3=1.0/float(ngrid*ngrid*ngrid)
      dtat=dt*rat
!
! move particles which fall in and belong to this processor : 
!		  use Lax-Wendroff second order two steps scheme
!                 1)  t^n ---> t^{n+1/2}
!                 2)  t^{n+1/2} ---> t^{n+1}
!
      do j=1,nparmax
!
	i3=int(x3(j)+0.5)	
	i3aux=i3
	if(i3.eq.0)i3=npes*nz
	ipe3=(i3-1)/nz
	if(ipe3.eq.mype)then
	  i3=i3-mype*nz
	  i1=int(x1(j)+0.5)
	  i2=int(x2(j)+0.5)
	  i1aux=i1
	  i2aux=i2
	  if(i1.eq.0)i1=nx
	  if(i2.eq.0)i2=ny
!
	  dd1=x1(j)-float(i1aux)+0.5
	  dd2=x2(j)-float(i2aux)+0.5
	  dd3=x3(j)-float(i3aux)+0.5
	  if(dd1.gt.1.0.or.dd1.lt.0.0)write(*,*)'WARNING: dd1 overflow'
	  if(dd2.gt.1.0.or.dd2.lt.0.0)write(*,*)'WARNING: dd2 overflow'
	  if(dd3.gt.1.0.or.dd3.lt.0.0)write(*,*)'WARNING: dd3 overflow'
	  de1=1.0-dd1
	  de2=1.0-dd2
	  de3=1.0-dd3
          d1= de1*de2*de3
          d2= dd1*de2*de3
          d3= de1*dd2*de3
          d4= dd1*dd2*de3
          d5= de1*de2*dd3
          d6= dd1*de2*dd3
          d7= de1*dd2*dd3
          d8= dd1*dd2*dd3
!
	  i1=i1
	  ip1=i1+1
	  if(i1.eq.nx)then
	    ip1=1
	  endif
!
	  j1=i2
	  jp1=i2+1
	  if(i2.eq.ny)then
	    jp1=1
	  endif
!
	  k1=i3
          kp1=i3+1
!
	  ffx=gxold(i1,j1,k1)*d1+&
              gxold(ip1,j1,k1)*d2+&
              gxold(i1,jp1,k1)*d3+&
              gxold(ip1,jp1,k1)*d4+&
              gxold(i1,j1,kp1)*d5+&
              gxold(ip1,j1,kp1)*d6+&
              gxold(i1,jp1,kp1)*d7+&
              gxold(ip1,jp1,kp1)*d8

	  ffy=gyold(i1,j1,k1)*d1+&
              gyold(ip1,j1,k1)*d2+&
              gyold(i1,jp1,k1)*d3+&
              gyold(ip1,jp1,k1)*d4+&
              gyold(i1,j1,kp1)*d5+&
              gyold(ip1,j1,kp1)*d6+&
              gyold(i1,jp1,kp1)*d7+&
              gyold(ip1,jp1,kp1)*d8

	  ffz=gzold(i1,j1,k1)*d1+&
              gzold(ip1,j1,k1)*d2+&
              gzold(i1,jp1,k1)*d3+&
              gzold(ip1,jp1,k1)*d4+&
              gzold(i1,j1,kp1)*d5+&
              gzold(ip1,j1,kp1)*d6+&
              gzold(i1,jp1,kp1)*d7+&
              gzold(ip1,jp1,kp1)*d8
!
! update velocity and position of the particle at t^n+1/2, and set periodic 
! boundary conditions
!
          vhalf(1)=v1(j)*(1.0-0.5*dat*dtat)+0.5*ffx*dtat
          vhalf(2)=v2(j)*(1.0-0.5*dat*dtat)+0.5*ffy*dtat
          vhalf(3)=v3(j)*(1.0-0.5*dat*dtat)+0.5*ffz*dtat
!
          xhalf(1)=x1(j)+0.5*v1(j)*dtat
          xhalf(2)=x2(j)+0.5*v2(j)*dtat
          xhalf(3)=x3(j)+0.5*v3(j)*dtat
	  do iidim1=1,2
            if (xhalf(iidim1).lt.0.0) xhalf(iidim1)=xhalf(iidim1)+xmax
            if (xhalf(iidim1).ge.xmax) xhalf(iidim1)=xhalf(iidim1)-xmax
	  enddo
!
	  i3=int(xhalf(3)+0.5)
	  i3aux=i3
	  if(i3.eq.0.and.mype.ne.0)i3=npes*nz
	  if(i3.eq.1.and.mype.ne.0)i3=npes*nz+1
	  i3=i3-mype*nz
          i1=int(xhalf(1)+0.5)
          i2=int(xhalf(2)+0.5)
	  i1aux=i1
	  i2aux=i2
          if(i1.eq.0)i1=nx
          if(i2.eq.0)i2=ny
!
          dd1=xhalf(1)-float(i1aux)+0.5
          dd2=xhalf(2)-float(i2aux)+0.5
          dd3=xhalf(3)-float(i3aux)+0.5
	  if(dd1.gt.1.0.or.dd1.lt.0.0)write(*,*)'WARNING: dd1 overflow'
	  if(dd2.gt.1.0.or.dd2.lt.0.0)write(*,*)'WARNING: dd2 overflow'
	  if(dd3.gt.1.0.or.dd3.lt.0.0)write(*,*)'WARNING: dd3 overflow'
          de1=1.0-dd1
          de2=1.0-dd2
          de3=1.0-dd3
          d1= de1*de2*de3
          d2= dd1*de2*de3
          d3= de1*dd2*de3
          d4= dd1*dd2*de3
          d5= de1*de2*dd3
          d6= dd1*de2*dd3
          d7= de1*dd2*dd3
          d8= dd1*dd2*dd3
!
          i1=i1
          ip1=i1+1
          if(i1.eq.nx)then
            ip1=1
          endif
!
          j1=i2
          jp1=i2+1
          if(i2.eq.ny)then
            jp1=1
          endif
!
          k1=i3
          kp1=i3+1
!
	  ffx=gx(i1,j1,k1)*d1+&
              gx(ip1,j1,k1)*d2+&
              gx(i1,jp1,k1)*d3+&
              gx(ip1,jp1,k1)*d4+&
              gx(i1,j1,kp1)*d5+&
              gx(ip1,j1,kp1)*d6+&
              gx(i1,jp1,kp1)*d7+&
              gx(ip1,jp1,kp1)*d8

          ffy=gy(i1,j1,k1)*d1+&
              gy(ip1,j1,k1)*d2+&
              gy(i1,jp1,k1)*d3+&
              gy(ip1,jp1,k1)*d4+&
              gy(i1,j1,kp1)*d5+&
              gy(ip1,j1,kp1)*d6+&
              gy(i1,jp1,kp1)*d7+&
              gy(ip1,jp1,kp1)*d8

          ffz=gz(i1,j1,k1)*d1+&
              gz(ip1,j1,k1)*d2+&
              gz(i1,jp1,k1)*d3+&
              gz(ip1,jp1,k1)*d4+&
              gz(i1,j1,kp1)*d5+&
              gz(ip1,j1,kp1)*d6+&
              gz(i1,jp1,kp1)*d7+&
              gz(ip1,jp1,kp1)*d8
!
! update velocity and position of the particle at t^{n+1}, and set periodic 
! boundary conditions
!
          v1(j)=v1(j)-rdtath*dath*vhalf(1)+rdtath*ffx
          v2(j)=v2(j)-rdtath*dath*vhalf(2)+rdtath*ffy
          v3(j)=v3(j)-rdtath*dath*vhalf(3)+rdtath*ffz
!
          x1(j)=x1(j)+rdtath*vhalf(1)
          x2(j)=x2(j)+rdtath*vhalf(2)
          x3(j)=x3(j)+rdtath*vhalf(3)
          if (x1(j).lt.0.0) x1(j)=x1(j)+xmax
          if (x1(j).ge.xmax) x1(j)=x1(j)-xmax
          if (x2(j).lt.0.0) x2(j)=x2(j)+xmax
          if (x2(j).ge.xmax) x2(j)=x2(j)-xmax
          if (x3(j).lt.0.0) x3(j)=x3(j)+xmax
          if (x3(j).ge.xmax) x3(j)=x3(j)-xmax
!
        endif
!
      enddo
!
! move particles that belong to other processors but fall in this one
!
      do j=1,n_tot
        i3=int(x3_aux_r(j)+0.5)
	i3aux=i3
	if(i3.eq.0)i3=nz*npes
	i3=i3-mype*nz
	if(i3.lt.0)i3=i3+(mype+1)*nz
        i1=int(x1_aux_r(j)+0.5)
        i2=int(x2_aux_r(j)+0.5)
	i1aux=i1
	i2aux=i2
        if(i1.eq.0)i1=nx
        if(i2.eq.0)i2=ny
!
        dd1=x1_aux_r(j)-float(i1aux)+0.5
        dd2=x2_aux_r(j)-float(i2aux)+0.5
        dd3=x3_aux_r(j)-float(i3aux)+0.5
	if(dd1.gt.1.0.or.dd1.lt.0.0)write(*,*)'WARNING: dd1 overflow'
	if(dd2.gt.1.0.or.dd2.lt.0.0)write(*,*)'WARNING: dd2 overflow'
	if(dd3.gt.1.0.or.dd3.lt.0.0)write(*,*)'WARNING: dd3 overflow'
        de1=1.0-dd1
        de2=1.0-dd2
        de3=1.0-dd3
        d1= de1*de2*de3
        d2= dd1*de2*de3
        d3= de1*dd2*de3
        d4= dd1*dd2*de3
        d5= de1*de2*dd3
        d6= dd1*de2*dd3
        d7= de1*dd2*dd3
        d8= dd1*dd2*dd3
!
        i1=i1
        ip1=i1+1
        if(i1.eq.nx)then
          ip1=1
        endif
!
        j1=i2
        jp1=i2+1
        if(i2.eq.ny)then
          jp1=1
        endif
!
        k1=i3
        kp1=i3+1
!
          ffx=gxold(i1,j1,k1)*d1+&
              gxold(ip1,j1,k1)*d2+&
              gxold(i1,jp1,k1)*d3+&
              gxold(ip1,jp1,k1)*d4+&
              gxold(i1,j1,kp1)*d5+&
              gxold(ip1,j1,kp1)*d6+&
              gxold(i1,jp1,kp1)*d7+&
              gxold(ip1,jp1,kp1)*d8

          ffy=gyold(i1,j1,k1)*d1+&
              gyold(ip1,j1,k1)*d2+&
              gyold(i1,jp1,k1)*d3+&
              gyold(ip1,jp1,k1)*d4+&
              gyold(i1,j1,kp1)*d5+&
              gyold(ip1,j1,kp1)*d6+&
              gyold(i1,jp1,kp1)*d7+&
              gyold(ip1,jp1,kp1)*d8

          ffz=gzold(i1,j1,k1)*d1+&
              gzold(ip1,j1,k1)*d2+&
              gzold(i1,jp1,k1)*d3+&
              gzold(ip1,jp1,k1)*d4+&
              gzold(i1,j1,kp1)*d5+&
              gzold(ip1,j1,kp1)*d6+&
              gzold(i1,jp1,kp1)*d7+&
              gzold(ip1,jp1,kp1)*d8
!
! update velocity and position of the particle at t^n+1/2, and set periodic
! boundary conditions
!
          vhalf(1)=v1_aux_r(j)*(1.0-0.5*dat*dtat)+0.5*ffx*dtat
          vhalf(2)=v2_aux_r(j)*(1.0-0.5*dat*dtat)+0.5*ffy*dtat
          vhalf(3)=v3_aux_r(j)*(1.0-0.5*dat*dtat)+0.5*ffz*dtat
!
          xhalf(1)=x1_aux_r(j)+0.5*v1_aux_r(j)*dtat
          xhalf(2)=x2_aux_r(j)+0.5*v2_aux_r(j)*dtat
          xhalf(3)=x3_aux_r(j)+0.5*v3_aux_r(j)*dtat
          do iidim1=1,2
            if (xhalf(iidim1).lt.0.0) xhalf(iidim1)=xhalf(iidim1)+xmax
            if (xhalf(iidim1).ge.xmax) xhalf(iidim1)=xhalf(iidim1)-xmax
          enddo
!
          i3=int(xhalf(3)+0.5)
	  i3aux=i3
	  if(mype.ne.0)then
	    if(i3.eq.0)i3=nz*npes
            i3=i3-mype*nz
	  endif
	  if(i3.lt.0)i3=i3+(mype+1)*nz
          i1=int(xhalf(1)+0.5)
          i2=int(xhalf(2)+0.5)
	  i1aux=i1
	  i2aux=i2
          if(i1.eq.0)i1=nx
          if(i2.eq.0)i2=ny
!
          dd1=xhalf(1)-float(i1aux)+0.5
          dd2=xhalf(2)-float(i2aux)+0.5
          dd3=xhalf(3)-float(i3aux)+0.5
	  if(dd1.gt.1.0.or.dd1.lt.0.0)write(*,*)'WARNING: dd1 overflow'
	  if(dd2.gt.1.0.or.dd2.lt.0.0)write(*,*)'WARNING: dd2 overflow'
	  if(dd3.gt.1.0.or.dd3.lt.0.0)write(*,*)'WARNING: dd3 overflow'
          de1=1.0-dd1
          de2=1.0-dd2
          de3=1.0-dd3
          d1= de1*de2*de3
          d2= dd1*de2*de3
          d3= de1*dd2*de3
          d4= dd1*dd2*de3
          d5= de1*de2*dd3
          d6= dd1*de2*dd3
          d7= de1*dd2*dd3
          d8= dd1*dd2*dd3
!
          i1=i1
          ip1=i1+1
          if(i1.eq.nx)then
            ip1=1
          endif
!
          j1=i2
          jp1=i2+1
          if(i2.eq.ny)then
            jp1=1
          endif
!
          k1=i3
          kp1=i3+1
!
          ffx=gx(i1,j1,k1)*d1+&
              gx(ip1,j1,k1)*d2+&
              gx(i1,jp1,k1)*d3+&
              gx(ip1,jp1,k1)*d4+&
              gx(i1,j1,kp1)*d5+&
              gx(ip1,j1,kp1)*d6+&
              gx(i1,jp1,kp1)*d7+&
              gx(ip1,jp1,kp1)*d8
          ffy=gy(i1,j1,k1)*d1+&
              gy(ip1,j1,k1)*d2+&
              gy(i1,jp1,k1)*d3+&
              gy(ip1,jp1,k1)*d4+&
              gy(i1,j1,kp1)*d5+&
              gy(ip1,j1,kp1)*d6+&
              gy(i1,jp1,kp1)*d7+&
              gy(ip1,jp1,kp1)*d8
          ffz=gz(i1,j1,k1)*d1+&
              gz(ip1,j1,k1)*d2+&
              gz(i1,jp1,k1)*d3+&
              gz(ip1,jp1,k1)*d4+&
              gz(i1,j1,kp1)*d5+&
              gz(ip1,j1,kp1)*d6+&
              gz(i1,jp1,kp1)*d7+&
              gz(ip1,jp1,kp1)*d8
! 
! update velocity and position of the particle at t^{n+1}, and set periodic
! boundary conditions
!
          v1_aux_r(j)=v1_aux_r(j)-rdtath*dath*vhalf(1)+rdtath*ffx
          v2_aux_r(j)=v2_aux_r(j)-rdtath*dath*vhalf(2)+rdtath*ffy
          v3_aux_r(j)=v3_aux_r(j)-rdtath*dath*vhalf(3)+rdtath*ffz
!
          x1_aux_r(j)=x1_aux_r(j)+rdtath*vhalf(1)
          x2_aux_r(j)=x2_aux_r(j)+rdtath*vhalf(2)
          x3_aux_r(j)=x3_aux_r(j)+rdtath*vhalf(3)
          if (x1_aux_r(j).lt.0.0) x1_aux_r(j)=x1_aux_r(j)+xmax
          if (x1_aux_r(j).ge.xmax) x1_aux_r(j)=x1_aux_r(j)-xmax
          if (x2_aux_r(j).lt.0.0) x2_aux_r(j)=x2_aux_r(j)+xmax
          if (x2_aux_r(j).ge.xmax) x2_aux_r(j)=x2_aux_r(j)-xmax
          if (x3_aux_r(j).lt.0.0) x3_aux_r(j)=x3_aux_r(j)+xmax
          if (x3_aux_r(j).ge.xmax) x3_aux_r(j)=x3_aux_r(j)-xmax
!
	enddo
!
END SUBROUTINE nbody

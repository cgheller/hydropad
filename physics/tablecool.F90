!
SUBROUTINE tablecool
!
        USE dimension
        USE matrix
        USE vector
        USE scalar
        USE ccool
        USE mpi_inc
        IMPLICIT NONE
        
! *******************
! * local variables *
! *******************

        integer :: i,j,k

! controllo del time step di sub-integrazione (cfr box sopra)
        real*8 :: dtfact, maxK, dtinf
        parameter (maxK=0.01, dtfact=0.25, dtinf=0.0)

! tempo cumulativo impiegato nell integrazione
        real*8 :: dttot

! varibili machine-dependent per i limiti di over/under flow
        real*8 :: maxd, mind, csi

! limiti machine-dependent per evitare over/under flow
        real*8 :: minlambda, maxlambda

! cooling time, cooling integration adaptive time-step
        real*8 :: tcool, dtcool


        real*8 :: inlambda,maxlambdat,minlambdat,meanlambdat
        real*8 :: maxflambdat,minflambdat,meanflambdat


! 
! CONTROL VARIABLES
!
        REAL*8  :: ncell, tcoolmin, tcoolmean, nsubmean, Deltat, Dp
        REAL*8  :: robustc, robc0, nsubmax 
        INTEGER :: steps

        REAL*8, DIMENSION(13) :: nsubsteps
        REAL*8, DIMENSION(19) :: dtth, dpph
           
        REAL*8  :: ncell_w, tcoolmin_w, tcoolmean_w, nsubmean_w
        REAL*8  :: robustc_w, robc0_w, nsubmax_w 

        REAL*8, DIMENSION(13) :: nsubsteps_w
        REAL*8, DIMENSION(19) :: dtth_w, dpph_w

        INTEGER :: integer_a, integer_b

!

        real*8 :: rhotot,rhoaux, tttaux, eaux, pold

        real*8 :: lambdat,Gettab       
        real*8 :: ph, lambdath, temph
       
! old variables
        
        integer :: downtable,uptable,downctable,upctable

        common/outtables/downtable,uptable,downctable,upctable

!     **********************
!     * cooling main cycle *
!     **********************

	dtth=0.0 
	dpph=0.0
	nsubsteps=0.0
	nsubmean=0.0
        maxd= huge(maxd)
        mind= tiny(mind)
        csi= maxd*mind

        ncell= 0
        robustc= 0.0
        robc0= 0.0
        tcoolmin= huge(maxd)
	tcoolmean=0.0
        nsubmax= 0.0


        downtable=0
        uptable=0
        downctable=0
        upctable=0
! @@@ outer do


        do k=1,nz
           do j=1,ny
              do i=1,nx

                ncell=ncell+1
                steps=0         

                 tttaux= tfact*p3d(i,j,k)/rho3d(i,j,k)
                 Deltat= tttaux
                 Dp= p3d(i,j,k)
                 rhoaux= rho3d(i,j,k)

                 lambdat= Gettab(rhoaux,tttaux)
                 inlambda=lambdat

                 dttot=0.0
                 Dp= p3d(i,j,k)
                 Deltat= ttt(i,j,k)
                 
! @@ inner do  ********

                 do
                    if(lambdat.eq.0)exit
                    if(dttot.ge.dt) exit
                
                    steps=steps+1

                    eaux= p3d(i,j,k)/(gamma-1.0)


                    if(eaux.lt.csi) then
                       minlambda= 0
                       maxlambda=maxd
                    else
                       minlambda= 1.0/(maxd/(eaux+eaux/10.0))
                       maxlambda= maxd/(eaux+eaux/10.0)
                    endif


                    if((lambdat.le.minlambda).and.&
                       (lambdat.ge.0)) then
                       tcool= maxd
                       dtcool= dtfact*dt
		

                    elseif((lambdat.gt.minlambda).and.&
                           (lambdat.lt.maxlambda))then
                       tcool= eaux/lambdat


!>>>>>>>>>>>>>>
!CONTROL::  store minimum and mean tcool
                                if(tcoolmin.gt.tcool) tcoolmin=tcool
                                tcoolmean=(tcoolmean*(real(ncell-1))+tcool)/real(ncell)
!<<<<<<<<<<<<<<

                       dtcool= min(maxK*tcool, dtfact*dt)
                       if(dtcool.lt.(dtinf*dt))dtcool= 0

                    elseif(lambdat.gt.maxlambda)then
                       tcool=0
                       dtcool=0
! >>>>>>>>>>
                       robc0=robc0+1.0
! <<<<<<<<<<
                    else
                       write(*,*)' --> help me at ',i,j,k,&
                                 ' : negative emissivity !!',mype 
                    endif
                    
                    if(dtcool.gt.(dt-dttot)) dtcool=(dt-dttot)
               

                    if(dtcool.gt.0) then

! ************* Lax-Wendroff

! i+1/2 guess
                      ph= p3d(i,j,k)-lambdat*dtcool*0.5*(gamma-1.0)
                      temph= tfact*ph/rho3d(i,j,k)
                      lambdath= Gettab(rhoaux,temph)

! i+1 
                      pold= p3d(i,j,k)
                      p3d(i,j,k)=p3d(i,j,k)-lambdath*dtcool*(gamma-1.0)
                      tttaux=ttt(i,j,k)
                      ttt(i,j,k)=tfact*p3d(i,j,k)/rho3d(i,j,k)
                    endif

                    if((dtcool.eq.0).or.(p3d(i,j,k).le.0)) then

                    if(p3d(i,j,k).le.0) then
                       write(*,*)mype
                       write(*,'(A,3(1x,I3),A)')& 
                       '>> unexpected robust cooling at ',i,j,k,&
                       ' p<=0  -> T down to ',k2ev,'K (1eV)'
                       write(*,*)'   physics at the cell is '
                       write(*,'(A,3(1x,E14.8))')'   p,t,r : ',&
                                    pold,tttaux,rho3d(i,j,k)
                       write(*,'(A,2(1x,E14.8))')'   lambda,lambdah: ',&
                                    lambdat,lambdath
                       write(*,'(A,2(1x,E14.8))')'   tcool,dtcool: ',&
                                    tcool,dtcool
          
                    else
! >>>>>>>>>>>> 
                       robustc=robustc+1.0
! <<<<<<<<<<<<
                    endif

                    ttt(i,j,k)= K2eV
                    p3d(i,j,k)=(ttt(i,j,k)/tfact)*rho3d(i,j,k)
                    dttot=dt
                    endif


                    rhoaux=rho3d(i,j,k)
                    tttaux=ttt(i,j,k)
                     
                    dttot=dttot+dtcool
                    if(dttot.ge.dt) exit


                    lambdat= Gettab(rhoaux,tttaux)
		 
! @@ end inner do  ********
                    
                 enddo

! >>>>>>>>>>>>>>>>
! CONTROL ::
                if(nsubmax.lt.real(steps)) nsubmax=real(steps)
                nsubmean= (nsubmean*(real(ncell-1))+real(steps))/real(ncell)

                 if(Deltat.gt.0) Deltat= (Deltat-ttt(i,j,k))/Deltat
                 
                 integer_a= int(Deltat*10)
                 if(integer_a.ge.1) then
                        dtth(integer_a)=dtth(integer_a)+1.0
                 else
                        integer_b= int(Deltat*100)-integer_a*10
                        dtth(10+integer_b)= dtth(10+integer_b)+1.0
                 endif
                 
                 if(Dp.gt.0) Dp= (Dp-p3d(i,j,k))/Dp
                 integer_a= int(Dp*10)
                 if(integer_a.ge.1) then
                        dpph(integer_a)=dpph(integer_a)+1.0
                 else
                        integer_b= int(Dp*100)-integer_a*10
                        dpph(10+integer_b)= dpph(10+integer_b)+1.0
                 endif

                 integer_a= steps/10
                 if(integer_a.lt.12) then
                        nsubsteps(integer_a+1)= nsubsteps(integer_a+1)+1.0
                 else
                        nsubsteps(13)= nsubsteps(13)+1.0
                 endif

! <<<<<<<<<<<<<<<<
! END CONTROL
                 
              enddo
           enddo
        enddo



        CALL MPI_Reduce(tcoolmin,tcoolmin_w,1,MPI_DOUBLE_PRECISION,&
                        MPI_MIN,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(tcoolmean,tcoolmean_w,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(nsubmean,nsubmean_w,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(nsubmax,nsubmax_w,1,MPI_DOUBLE_PRECISION,&
                        MPI_MAX,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(robustc,robustc_w,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(robc0,robc0_w,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(dtth,dtth_w,19,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(dpph,dpph_w,19,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Reduce(nsubsteps,nsubsteps_w,13,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(mype.eq.0)then
          tcoolmean_w=tcoolmean_w/real(npes)
          nsubmean_w=nsubmean_w/real(npes)
          write(600,200)redshift,tcoolmin_w,tcoolmean_w,dt,&
	                nsubmean_w,nsubmax_w,robustc_w,robc0_w
          write(601,*)redshift
          do i=1,19
            write(601,201)i,dtth_w(i),dpph_w(i)
          enddo
          do i=1,13
            write(601,202)i,nsubsteps_w(i)
          enddo
200       format(8(1x,e13.7))
201       format(1x,i2,2(1x,e13.7))
202       format(1x,i2,1x,e13.7)
        endif


        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE tablecool

!
!********************************************************************
!

      real*8 function Gettab(bmdens,ttemp)
        USE dimension
        USE vector
        USE scalar
        USE ccool
        IMPLICIT NONE

      integer GetIndexTable

      real*8 mr,qr,mr2,qr2,mr3,qr3
      real*8 bmdens,ttemp

      integer indexr,indexr2,indexr3,searchst
      integer downtable,uptable,downctable,upctable

      common/outtables/downtable,uptable,downctable,upctable

! cerca l intervallo di temperatura di appartenenza

      indexr=  GetIndexTable(1,ttemp,searchst)
      indexr2= GetIndexTable(2,ttemp,searchst)
      indexr3= GetIndexTable(3,ttemp,searchst)      


! calcola la retta per l intervallo individuato

      if(indexr.gt.0)then
         mr= mrtable(indexr); qr=qrtable(indexr)
      elseif(indexr.eq.-2)then
         mr= 0.0; qr=mrtable(rhotable)*trtable(rhotable)+&
                  qrtable(rhotable)
         uptable=uptable+1
      else
         mr=0; qr=0
         downtable=downtable+1
      endif

      if(indexr2.gt.0)then
         mr2= mr2table(indexr2); qr2=qr2table(indexr2)
      elseif(indexr2.eq.-2)then
         mr2= 0.0; qr2=mr2table(rho2table)*tr2table(rho2table)+&
                  qr2table(rho2table)
         upctable=upctable+1
      else
         mr2=0; qr2=0
         downctable=downctable+1
      endif

      if(indexr3.gt.0)then
         mr3= mr3table(indexr3); qr3=qr3table(indexr3)
      elseif(indexr3.eq.-2)then
         mr3= 0.0; qr3=mr3table(rho3table)*tr3table(rho3table)+&
                  qr3table(rho3table)
         upctable=upctable+1
      else
         mr3=0; qr3=0
         downctable=downctable+1
      endif

! calcola l emissione sommando i singoli contributi



      Gettab= (mr*ttemp +qr)*&
			 bmdens * cmfact/hfrac/at**4 +&
	      (mr2*ttemp +qr2)*&
                         bmdens**2 * hfact/hfrac**2/at**3 +&
              (mr3*ttemp +qr3)*&
                         bmdens**3 * hfact/hfrac**2 * &
                         rhonow / mh / at**6

 
END FUNCTION Gettab

! ******************
! * binary search  *
! ******************

! search tval in matrix ttable that have ntable elements 
! (matrix is ascending ordered)      

      integer function GetIndexTable(q,tval,searchst)

        USE dimension
        USE vector
        USE scalar
        USE ccool
        IMPLICIT NONE 

      real*8 tval,vh
      real*8 maxval,minval

      integer n,q,searchst
      integer found
      integer half,upi,downi
      integer result,greatest

      downi=1

      if(q.eq.1) then
! search in temp table of rho**1 cooling curve

         greatest=rhotable
         maxval= trtable(rhotable)
         minval= trtable(1)

      else if(q.eq.2) then
! search in temp table of rho**2 cooling curve

         greatest=rho2table
         maxval= tr2table(rho2table)
         minval= tr2table(1)

      else
! search in temp table of rho**3 cooling curve

         greatest=rho3table
         maxval= tr3table(rho3table)
         minval= tr3table(1)
 
      endif

      upi=greatest
      found=0
      searchst=0

      if(tval.gt.maxval) then
         result= -2
      elseif(tval.lt.minval)then
         result= -3
      else

      do 100
         if((upi-downi).eq.1) found=1
         if(found.eq.1)exit

         half=(upi+downi)/2
         if(q.eq.1)then
            vh=trtable(half)
         else if(q.eq.2)then
            vh=tr2table(half)
         else 
            vh=tr3table(half)
         endif

         if(tval.lt.vh) then
            upi=half
         elseif(tval.gt.vh) then
            downi=half
         else
            downi=half
            found=1
         endif
         searchst=searchst+1
 100  continue

      result=downi

      endif

      GetIndexTable= result

END FUNCTION GetIndexTable

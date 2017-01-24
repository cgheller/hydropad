FUNCTION GetFitPoints(ndata,xvector,vvector,fpoints)
      implicit None
      integer :: GetFitPoints
      integer :: ndata,cfp,i
      integer :: fpoints(ndata)
      real*8  :: xvector(ndata),vvector(ndata)
      
      integer :: xlast,xnow,half,xmaxdiff
      real*8  :: linfit
      real*8  :: m,diff,maxdiff,der1,der2
      real*8  :: maxv,minv,gap
      real*8  :: r1,r2,m1,m2,q1,q2
      real*8  :: x1(ndata),x2(ndata),y1(ndata),y2(ndata)
      integer :: jump
      real*8  :: maxp,frac
      parameter(jump=1,maxp=0.01,frac=0.02)
      
    cfp=1
      fpoints(1)=1
      xlast=1
      xnow=2

      do 11
         
         if(xnow.eq.ndata) exit

         m= (vvector(xnow)-vvector(xlast))/&
            (xvector(xnow)-xvector(xlast))
    
         half= (xnow+xlast)/2

         maxdiff=0
        
         
         if((xnow-xlast).gt.1) then
             maxv=-huge(maxv)
             minv=huge(maxv)
             do i=xlast+1,xlast+(xnow-xlast-1)
                maxv= max(maxv,vvector(i))
                minv= min(minv,vvector(i))
             enddo
                gap= abs(maxv-minv)
             if(gap.ne.0)then
                do i=xlast+1,xlast+(xnow-xlast-1)
                   diff= abs(&
                        (m*xvector(i)+(vvector(xlast)-m*xvector(xlast))- &
                        vvector(i))/ gap )
                   if(diff.gt.maxdiff) xmaxdiff=i
                   maxdiff=max(maxdiff,diff)
                enddo
             endif
          endif

         if(maxdiff.gt.maxp) then
            if(((float(xnow-xlast)/float(ndata)).gt.frac).and.&
            ((xnow-xlast).gt.4).and.((xnow-xmaxdiff).gt.2)) then
               if(xmaxdiff.eq.(xlast+1))xmaxdiff=xlast+2
               cfp=cfp+1
               fpoints(cfp)=xmaxdiff
               xlast=xmaxdiff
            endif
            if((xnow-xlast).gt.3)then
               xlast=xnow-2
            elseif((xnow-xlast).eq.3) then
               xlast=xnow-1
            endif
            xnow=xlast+1
            cfp=cfp+1
            fpoints(cfp)=xlast

         else
            xnow=xnow+jump
            if(xnow.gt.ndata)xnow=ndata
         endif
            
 11   continue

      if((xnow-xlast).lt.3) then
         fpoints(cfp)=ndata
      else
       cfp=cfp+1
       fpoints(cfp)=ndata
      endif

      GetFitPoints=cfp

END FUNCTION GetFitPoints


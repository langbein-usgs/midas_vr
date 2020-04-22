      PROGRAM MIDAS_VR
C  MODIFIED from MIDAS by John Langbein, Sept 2019.
c     the VR -variance in rate, added by Langbein (JL)
c
c
c   Author: Geoff Blewitt.  Copyright (C) 2015.
c
c   (1) take all possible pairs of data seperated by precisely one year.
c   (2) if a corresponding pair cannot be found, match with the next 
c       closest epoch that has not yet been used. 
c       2015-07-17: can have unused epochs; particularly bad for campaigns.
c       Solution: make algorithm symmetric, forward and back in time.
c       2015-10-10: MIDAS4 has option to use a step epoch file.  Pairs are
c       not allowed to cross epochs. Daily solution of step epoch is ignored.
c       Modify 
c   (3) find the median, v50, and median of absolute deviation (MAD)
c   (4) trim distribution for data > 2 std dev
c   (5) recompute median and MAD
c   (6) estimate standard error in median using MAD
C   (7) JL-add -- estimate standard error using the VR algorithm as #6 is not correct
c
c   Read data from MIDAS.TENV in tenv2 format, 
c    or from MIDAS.TENU in free format: station, time, east, north, up
c   Optionally read step epochs from MIDAS.STEPS
c   Write one line of statistics out to stdout, and write
c   to MIDAS.RENV a new tenv2 file containing residuals to the fit. 
c
c   JL -- Add subroutine to implement Hackl (2011) AVR routine to compute
c   rate errors; modified AVR to do medians, etc...
c
c
      implicit none

c   following is the label and version number for this fit
      character*9, parameter :: label = 'MIDAS_VR1'
      integer, parameter :: maxm=9999, minm=6, minn=10, maxn=19999
      integer, parameter :: maxstep=99 
c     setting tol=0.001d+0 forces pairs separated by precisely 1 year
      real*8,  parameter :: tol=0.001d+0, tolneg = -tol
      real*8,  parameter :: deltmin = 1.d0+real(minm)/365.25+tol
      real*8,  parameter :: svmax=9.999999d+0, rmax=99.999999d+0

      integer i, j, k, m, n1, n2, n, mgood, ndata, ios, nstep, istep
      logical freeformat
      integer ip(2,maxn) 
      real*8 delt, fdt, dt
      integer ne, nn, nu
      real*8 ve50, vn50, vu50, de50, dn50, du50, ce, cn, cu
      real*8 sve, svn, svu, sde, sdn, sdu
      real*8 fe, fn, fu
      real*8 fit_e,fit_n,fit_u
      real*8 ts, tstep(maxstep+1), tbstep(maxstep+1)
      real*8 t(maxm), tb(maxm)
      real*8 ve(maxn), vn(maxn), vu(maxn)
      real*8 de(maxn), dn(maxn), du(maxn)
      logical good(maxm)
      real*8 xe50, xn50, xu50 
      real Tstart,Tend
      character*80 line
      character*4 sta, stepsta
      character*7 date(maxm)
      real*8 xe(maxm), xn(maxm), xu(maxm)
      real*8 se(maxm), sn(maxm), su(maxm), ah(maxm)
      real*8 c1(maxm), c2(maxm), c3(maxm) 
      real*8 re(maxm), rn(maxm), ru(maxm)
      integer mjd(maxm), mgpsw(maxm), mwday(maxm)
      real*8 qmedian,twind
      real*8 rate_error_e,rate_error_n, rate_error_u

      Tstart=second()
c-----DATA INPUT SECTION---------------------------------------------

      open(unit=91,action='write',file='MIDAS.ERR',form='formatted',
     +     iostat=ios, status='new')
      open(unit=51,action='read',file='MIDAS.TENV',form='formatted',
     +     iostat=ios, status='old') 
      if(ios.eq.0) then
         freeformat = .FALSE.
      else
         open(unit=51,action='read',file='MIDAS.TENU',form='formatted',
     +     iostat=ios, status='old')
         if(ios.ne.0) then
           write(91,*) 'FATAL: Input file not found'
           call exit(1)
         endif
         freeformat = .TRUE.
      endif


      m = 1
      if(freeformat) then
        do
          read(51,*,iostat=ios) sta, t(m), xe(m), xn(m), xu(m)
          if(ios.ne.0) exit
          m = m + 1  
          if(m > maxm) then
             write(91,'(2a,i4)') sta,' WARNING: m epochs > ',maxm
             exit
          endif
        enddo
      else
        do
          read(51,*,iostat=ios) sta, date(m), t(m), mjd(m), mgpsw(m), 
     +      mwday(m), xe(m), xn(m), xu(m), ah(m), se(m), sn(m), su(m),
     +      c1(m), c2(m), c3(m)
          if(ios.ne.0) exit
          if(xe(m)>99.9) cycle
          if(xn(m)>99.9) cycle
          if(xu(m)>99.9) cycle
          if(se(m)>0.03) cycle
          if(sn(m)>0.03) cycle
          if(su(m)>0.1) cycle
          m = m + 1  
          if(m > maxm) then
             write(91,'(2a,i4)') sta,' WARNING: m epochs >= ',maxm
             exit
          endif
        enddo
      endif
      m = m - 1
      close(51)

      if ( m < 1 )  then
         write(91,'(2a)') sta,' FATAL: No data'
         call exit(1)
      endif

      if ( m < 2 )  then
         write(91,'(2a)') sta,' FATAL: No pairs'
         call exit(1)
      endif

c   Compute time span
      delt = t(m) - t(1)
      if ( delt < deltmin )  then
         write(91,'(2a,f5.3,a,f5.3)') sta,' FATAL: Time span(yr) ',
     +      delt,' < ',deltmin
         call exit(1)
      endif

      if ( m < minm ) then
         write(91,'(2a,i1)') sta,' FATAL: m epochs < ',minm
         call exit(1)
      endif

      nstep = 0
      open(unit=52,action='read',file='MIDAS.STEPIN',form='formatted',
     +     iostat=ios, status='old')
      if(ios.eq.0) then
         open(unit=62,action='write',file='MIDAS.STEPOUT',
     +      form='formatted',iostat=ios, status='new')
         do 
           read(52,'(a80)',iostat=ios) line
           if(ios.ne.0) exit
           read(line,*,iostat=ios) stepsta, ts
           if(ios.ne.0) exit

c          only use steps from this station within the data span
           if(stepsta.ne.sta) cycle
           if(ts<t(1) .or. ts>t(m) ) cycle

           if(nstep>=maxstep) then
             write(91,'(2a,i4)') sta,' WARNING: nstep > ',maxstep
             exit
           endif
           nstep = nstep + 1
           tstep(nstep) = ts
           write(62,'(a)') trim(line)
         enddo
      endif

c-----DATA SELECTION SECTION-----------------------------------------
      twind=1.0
c   select pairs forward in positive time
      call selectpair(m,maxn,tol,t,n1,ip,nstep,tstep,twind)
      
c   now select pairs backward in negative time to ensure 
c     the algorithm is time symmetric
      call tback(m,t,tb,nstep,tstep,tbstep)
      call selectpair(m,maxn,tol,tb,n2,ip(1,n1+1),nstep,tbstep,twind)
c
c     correct pair indices for forward positive time
      n = n1+n2
      do k = n1+1, n
         i = ip(1,k)
         j = ip(2,k)
         ip(1,k) = m+1-j
         ip(2,k) = m+1-i
      enddo

c   initialize data points contributing to median
      do i = 1, m
         good(i) = .FALSE.
      enddo

c   compute velocity data
      do k = 1, n
         i = ip(1,k)
         j = ip(2,k)
         good(i) = .TRUE.
         good(j) = .TRUE.
         dt = t(j) - t(i)

c         debug/scrutinize pair selection
c         print'(2i6,2a8,3f10.4)',i,j,date(i),date(j),t(i),t(j),dt

         ve(k) = (xe(j)-xe(i))/dt
         vn(k) = (xn(j)-xn(i))/dt
         vu(k) = (xu(j)-xu(i))/dt
      enddo

      if ( n >= maxn ) then
         write(91,'(2a,i5)') sta,' WARNING: n pairs >= ',maxn
      endif
      if ( n < minn )  then
         write(91,'(2a,i2)') sta,' FATAL: n pairs < ',minn
         call exit(1)
      endif

c   Number of data points used
      mgood = 0
      do i = 1, m
         if ( good(i) ) then
            mgood = mgood + 1
         endif
      enddo

c-----ESTIMATION SECTION---------------------------------------------

c   Median of the distribution
      ve50 = qmedian(n,ve)
      vn50 = qmedian(n,vn)
      vu50 = qmedian(n,vu)

c   Absolute velocity deviation from the median
      do i = 1, n
         de(i) = abs(ve(i)-ve50)
         dn(i) = abs(vn(i)-vn50)
         du(i) = abs(vu(i)-vu50)
      enddo
      
c   Median Absolute Deviation (MAD)
      de50 = qmedian(n,de)
      dn50 = qmedian(n,dn)
      du50 = qmedian(n,du)

c   Estimated standard deviation of velocities
      sde = 1.4826*de50
      sdn = 1.4826*dn50
      sdu = 1.4826*du50

c   Delete velocities more than 2 std deviations from median
      ce = 2.0*sde
      cn = 2.0*sdn
      cu = 2.0*sdu
      ne = 0
      nn = 0
      nu = 0
      do i = 1, n
         if (de(i) < ce) then
            ne = ne + 1
            ve(ne) = ve(i)
         endif
         if (dn(i) < cn) then
            nn = nn + 1
            vn(nn) = vn(i)
         endif
         if (du(i) < cu) then
            nu = nu + 1
            vu(nu) = vu(i)
         endif
      enddo

c     Recompute median
      ve50 = qmedian(ne,ve)
      vn50 = qmedian(nn,vn)
      vu50 = qmedian(nu,vu)

c     Recompute absolute deviations
      do i = 1, ne
         de(i) = abs(ve(i)-ve50)
      enddo
      do i = 1, nn
         dn(i) = abs(vn(i)-vn50)
      enddo
      do i = 1, nu
         du(i) = abs(vu(i)-vu50)
      enddo
      
c     Recompute MAD     
      de50 = qmedian(ne,de)
      dn50 = qmedian(nn,dn)
      du50 = qmedian(nu,du)

c     Estimated standard deviation of velocities
c     Multiply by theoretical factor of 1.4826 
      sde = 1.4826*de50
      sdn = 1.4826*dn50
      sdu = 1.4826*du50

c     Standard errors for the median velocity
c     Multiply by theoretical factor of sqrt(pi/2) = 1.2533
c     Divide number of data by 4 since we use coordinate data a nominal 4 times
      sve = 1.2533*sde/sqrt(real(ne)/4.0)
      svn = 1.2533*sdn/sqrt(real(nn)/4.0)
      svu = 1.2533*sdu/sqrt(real(nu)/4.0)

c     Scale standard errors by ad hoc factor of 3 to be realistic
      sve = 3.0*sve
      svn = 3.0*svn
      svu = 3.0*svu

c   Compute intercept at first epoch
c   Intercepts: xe50, xn50, xu50
      do i = 1, m
         dt = t(i) - t(1)
         re(i) = xe(i) - ve50*dt
         rn(i) = xn(i) - vn50*dt
         ru(i) = xu(i) - vu50*dt
      enddo
      xe50 = qmedian(m,re)
      xn50 = qmedian(m,rn)
      xu50 = qmedian(m,ru)
c
c   Compute residuals
      do i = 1, m
         re(i) = re(i) - xe50
         rn(i) = rn(i) - xn50
         ru(i) = ru(i) - xu50
      enddo

c-----OUTPUT RESULTS SECTION-----------------------------------------

c   Protect against format errors and identify gross problems
      if ( ve50 > svmax .or. ve50 < -svmax ) ve50 = svmax
      if ( vn50 > svmax .or. vn50 < -svmax ) vn50 = svmax
      if ( vu50 > svmax .or. vu50 < -svmax ) vu50 = svmax
      if ( sve > svmax ) sve = svmax
      if ( svn > svmax ) svn = svmax
      if ( svu > svmax ) svu = svmax
      if ( sde > svmax ) sde = svmax
      if ( sdn > svmax ) sdn = svmax
      if ( sdu > svmax ) sdu = svmax
      do i = 1, m
        if ( re(i) > rmax .or. re(i) < -rmax ) re(i) = rmax
        if ( rn(i) > rmax .or. rn(i) < -rmax ) rn(i) = rmax
        if ( ru(i) > rmax .or. ru(i) < -rmax ) ru(i) = rmax
      enddo

c   Fraction of pairs removed
      fe = real(n-ne)/real(n)
      fn = real(n-nn)/real(n)
      fu = real(n-nu)/real(n)
      Tend=second()
c      print*,' Midas cpu time',Tend-Tstart
      Tstart=Tend
c  Compute rate uncertainty based upon an approximation of Power law noise
      open(52,file="avr_e.dat")
      open(53,file="avr_e.fit")
      open(54,file="avr_e.pred")
      call AVR(t,re,m,maxm,nstep,tstep,rate_error_e,fit_e)
      close(52)
      close(53)
      close(54)
      open(52,file="avr_n.dat")
      open(53,file="avr_n.fit")
      open(54,file="avr_n.pred")
      call AVR(t,rn,m,maxm,nstep,tstep,rate_error_n,fit_n)
      close(52)
      close(53)
      close(54)
      open(52,file="avr_u.dat")
      open(53,file="avr_u.fit")
      open(54,file="avr_u.pred")
      call AVR(t,ru,m,maxm,nstep,tstep,rate_error_u,fit_u)
      close(52)
      close(53)
      close(54)

      Tend=second()
c      print*,' Midas VR cpu time',Tend-Tstart
c   Write out solution  
      write(6,'(a4,1x,a10,2f10.4,f8.4,i5,i5,i7,
     +        3f10.6,1x,3f9.6,
     +        3f11.6, 3f6.3, 3f9.6, i3,1x,3f9.6,1x,3(1x,f13.2))')
     +        sta, label, t(1), t(m), delt, m, mgood, n,
     +        ve50, vn50, vu50, sve, svn, svu,
     +        xe50, xn50, xu50, fe, fn, fu, sde, sdn, sdu, nstep,
     +        rate_error_e,rate_error_n,rate_error_u,
     +        fit_e,fit_n,fit_u
c
c   Write residuals
      if(freeformat) then
        open(unit=61,action='write',file='MIDAS.RENU',form='formatted',
     +     iostat=ios, status='new')
        if(ios.ne.0) then
           write(91,'(2a)') sta,' FATAL: Cannot open output file'
           call exit(1) 
        endif
        do i = 1, m
          write(61,'(a4,x,f9.4,3f11.6)')
     +       sta, t(i), re(i), rn(i), ru(i)
        enddo
      else
        open(unit=61,action='write',file='MIDAS.RENV',form='formatted',
     +     iostat=ios, status='new')
        if(ios.ne.0) then
           write(91,'(2a)') sta,' FATAL: Cannot open output file'
           call exit(1) 
        endif
        do i = 1, m
          write(61,'(a4,x,a7,x,f9.4,i6,i5,i2,3f11.6,f8.4,3f9.6,3f10.6)')
     +       sta, date(i), t(i), mjd(i), mgpsw(i), mwday(i), 
     +       re(i), rn(i), ru(i), ah(i), se(i), sn(i), su(i), 
     +       c1(i), c2(i), c3(i)
        enddo
      endif
      close(61)
c
c  call subroutine that computes better estimates of rate error
c


      end
c
c
      subroutine AVR(time,res,nobs,maxm,nstep,tstep,rate_error,fit_x)
      integer iday(maxm),iiday(maxm),wn_num,nobs,icol(7),istep(100)
      integer icount(1000)
      real*8 time(maxm),res(maxm),day0,sumavs,rres(maxm),sumvar
      real*8 tims(1000),rate(100000),chi2(1000),sigavr(1000),dum(1000)
      real*8 qmedian,wn,plamp1,plep,rsum,A(1000,7),tstep(100)
      real*8 sig_scale plexpFin wnFin rmmin,fit_x
      real*8 minav, maxav, rate_error,tlen
      real*8 x(5),e(5),periodic,xbest(5)
      character*5 type
      character*5 tbest
      character*1 MADy, AVRy, FLRWy,OUTy
c  MADy and AVRy are, at least for new, set in the compiled version
c   MADy  controls whether rate-variance is computed using MAD, or sumvar/n
c   AVRy  controls whether Allan Variance in rate (Hackl) or Variance in rate (more or less a jacknife)
c      MADy='y'
c      AVRy='y'
       open(66,file="avr.config")
       rewind(66)
       read(66,*)MADy
       read(66,*)AVRy
c  Xscale is scaling for apprior weighting of chi^2 for fitting.
       read(66,*)Xscale
c  Output stuff for plotting
       read(66,*)OUTy
       close (66)
c  Thres  is the threshold sigma to toss-out estimates of AVR or VR
       Thres=2.
       FLRWy='y'
c  compute day number from time in years
      do i=1,nobs
        if (i .eq. 1) day0=365.25*time(1)
        iday(i)=int(time(i)*365.25-day0+0.05+1.)
c        write(51,*)i,iday(i),time(i),res(i)
      end do
c  compute day number for step
      do i=1,nstep
        istep(i)=int(tstep(i)*365.25-day0+0.05+1.)
      end do
      stt = iday(1)
      ent = iday(nobs)
      minav = 7.d0
      maxav = (ent-stt)/4.d0
c  generate values of times to sample AVR that are equally spaced in
c     log(time)
      numav=int(((log(maxav) - log(minav))/0.075) +0.5 )
     &   + 1
      kk=0
      do k=1,numav
        timxs=int(exp(log(7.) + 0.075*(k-1)) +0.5)
        if (k .eq. 1) then
          timxsLast=timxs
          kk=kk+1
          tims(kk)=timxs
        else
          if (timxs .ne. timxsLast) then
            kk=kk+1
            timxsLast=timxs
            tims(kk)=timxs
          end if
        end if
      end do
      tims(kk+1)=maxav
      numav=kk+1
c      do k=1,numav
c      write(52,*)k,tims(k)
c      end do
c
c --------------------------------------------------------------------------
c
c  get so-called chi2 values
c
      if ((AVRy .eq. 'y') .and. (MADy .eq. 'y')) then
c
c  compute AVR for each bin with length tims
c
      do i=1,numav
        dt=tims(i)
c        print*,i,dt
        num = 0
        sumavs = 0.d0
        avr_last=1.0e+20
        do it = 0, nint((ent-stt)/dt)
          tx = stt + it*dt
          call getdata(tx,tx+dt,iday,res,iiday,rres,numdt,nobs,maxm,
     &      nstep,istep)
c          write(40,*)it,tx,tx+dt,numdt
          if (numdt .gt. (0.75*dt)) then
            sumy=0.
            sumt=0.
            do k=1,numdt
              sumt=sumt+iiday(k)
              sumy=sumy+rres(k)
            end do
            avet=sumt/float(numdt)
            avey=sumy/float(numdt)
            sumyt=0.0
            sumtt=0.0
            do k=1,numdt
              sumtt=sumtt+(iiday(k)-avet)**2
              sumyt=sumyt+(rres(k)-avey)*
     &                  (iiday(k)-avet)
            end do
            ratex=sumyt/sumtt
            if (avr_last .gt. .9e+20) then
               avr_last=ratex
            else 
               rate(num+1)=((ratex-avr_last)**2)/2.
               avr_last=ratex
               num=num+1
             end if
          end if
        end do
c  compute MAD of rate (note that rate is in absolute value here
c        chi2(i)=(sumavs/float(num))
        if (num .gt. 2) then
          xrate=qmedian(num,rate)
          yrate=xrate
          chi2(i)=(1.4826*xrate)*1.4826
c  remove rates that exceed 2 standard deviations
          rateLimit=Thres*xrate*(1.4826**2)
          kk=0
          do k=1,num
            if (rate(k) .le. rateLimit) then
             kk=kk+1
             rate(kk)=rate(k)
            end if
          end do
            num=kk
          if (num .gt. 2) then 
            xrate=qmedian(num,rate)
            chi2(i)=(1.4826*xrate)*1.4826  
c   note that results of fitting chi2 is dependent on selecting
c     the exponent of num^index ---  1 seems to work best....
            sigavr(i)=chi2(i)/(float(num)**Xscale)
            icount(i)=num
          else
            chi2(i)=1.0e+20
          end if
        else
          chi2(i)=1.0e+20
        end if
      end do
      end if
c  end ((AVRy .eq. 'y') .and. (MADy .eq. 'y'))

      if ((AVRy .eq. 'n') .and. (MADy .eq. 'y')) then
c
c  compute AVR for each bin with length tims
c
      nloop=2
c      print*,'nloop', nloop
      do i=1,numav
        dt=tims(i)
c        print*,i,dt
        num = 0
        sumavs = 0.d0
        do it = 0, nint((ent-stt)/dt)
        do iloop=1,nloop
          tx = stt + it*dt+(iloop-1)*int(dt/nloop)
          call getdata(tx,tx+dt,iday,res,iiday,rres,numdt,nobs,maxm,
     &      nstep,istep)
c          write(40,*)it,tx,tx+dt,numdt
          if (numdt .gt. (0.75*dt)) then
            sumy=0.
            sumt=0.
            do k=1,numdt
              sumt=sumt+iiday(k)
              sumy=sumy+rres(k)
            end do
            avet=sumt/float(numdt)
            avey=sumy/float(numdt)
            sumyt=0.0
            sumtt=0.0
            do k=1,numdt
              sumtt=sumtt+(iiday(k)-avet)**2
              sumyt=sumyt+(rres(k)-avey)*
     &                  (iiday(k)-avet)
            end do
            ratex=sumyt/sumtt
            rate(num+1)=ratex**2
            num=num+1
          end if
        end do
        end do
c  compute MAD of rate (note that rate is in absolute value here
c        chi2(i)=(sumavs/float(num))
        if (num .gt. 2) then
          xrate=qmedian(num,rate)
          yrate=xrate
          chi2(i)=(1.4826*xrate)*1.4826
c  remove rates that exceed 2 standard deviations
          rateLimit=Thres*xrate*(1.4826**2)
          kk=0
          do k=1,num
            if (rate(k) .le. rateLimit) then
             kk=kk+1
             rate(kk)=rate(k)
            end if
          end do
            num=kk
          if (num .gt. 2) then 
            xrate=qmedian(num,rate)
            chi2(i)=(1.4826*xrate)*1.4826  
c   note that results of fitting chi2 is dependent on selecting
c     the exponent of num^index ---  1 seems to work best....
            sigavr(i)=chi2(i)/(float(num)**Xscale)
            icount(i)=num
          else
            chi2(i)=1.0e+20
          end if
        else
          chi2(i)=1.0e+20
        end if
      end do
      end if
c  end ((AVRy .eq. 'n') .and. (MADy .eq. 'y')) 

      if ((AVRy .eq. 'y') .and. (MADy .eq. 'n')) then
c
c  compute AVR for each bin with length tims
c
      do i=1,numav
        dt=tims(i)
c        print*,i,dt
        num = 0
        sumavs = 0.d0
        avr_last=1.0e+20
        do it = 0, nint((ent-stt)/dt)
          tx = stt + it*dt
          call getdata(tx,tx+dt,iday,res,iiday,rres,numdt,nobs,maxm,
     &      nstep,istep)
c          write(40,*)it,tx,tx+dt,numdt
          if (numdt .gt. (0.75*dt)) then
            sumy=0.
            sumt=0.
            do k=1,numdt
              sumt=sumt+iiday(k)
              sumy=sumy+rres(k)
            end do
            avet=sumt/float(numdt)
            avey=sumy/float(numdt)
            sumyt=0.0
            sumtt=0.0
            do k=1,numdt
              sumtt=sumtt+(iiday(k)-avet)**2
              sumyt=sumyt+(rres(k)-avey)*
     &                  (iiday(k)-avet)
            end do
            ratex=sumyt/sumtt
            if (avr_last .gt. .9e+20) then
               avr_last=ratex
            else 
               rate(num+1)=((ratex-avr_last)**2)/2.
               sumavs=sumavs+rate(num+1)
               avr_last=ratex
               num=num+1
             end if
          end if
        end do
c        sumvar=0.0
c        do k=1,num
c          sumvar=sumvar+rate(k)
c        end do
c  compute MAD of rate (note that rate is in absolute value here
c        chi2(i)=(sumavs/float(num))
        if (num .gt. 2) then
          
          chi2(i)=sumavs/float(num)
c          print*,i,num,chi2(i)
c  remove rates that exceed 2 standard deviations
          rateLimit=Thres*chi2(i)*(1.4826**2)
          kk=0
          do k=1,num
            if (rate(k) .le. rateLimit) then
             kk=kk+1
             rate(kk)=rate(k)
            end if
          end do
            num=kk
          if (num .gt. 2) then 
            sumavs=0.0
            do kk=1,num
              sumavs=sumavs+rate(kk)
            end do

            chi2(i)=sumavs/float(num)
c            print*,'  ',i,num,chi2(i)

c   note that results of fitting chi2 is dependent on selecting
c     the exponent of num^index ---  1 seems to work best....
            sigavr(i)=chi2(i)/(float(num)**Xscale)
            icount(i)=num
          else
            chi2(i)=1.0e+20
          end if
        else
          chi2(i)=1.0e+20
        end if
      end do
      end if
c  end ((AVRy .eq. 'y') .and. (MADy .eq. 'n'))

      if ((AVRy .eq. 'n') .and. (MADy .eq. 'n')) then
c
c  compute AVR for each bin with length tims
c
      do i=1,numav
        dt=tims(i)
c        print*,i,dt
        num = 0
        sumavs = 0.d0
        nloop=2
        do it = 0, nint((ent-stt)/dt)
        do iloop=1,nloop
          tx = stt + it*dt+(iloop-1)*int(dt/nloop)
          call getdata(tx,tx+dt,iday,res,iiday,rres,numdt,nobs,maxm,
     &      nstep,istep)
c          write(40,*)it,tx,tx+dt,numdt
          if (numdt .gt. (0.75*dt)) then
            sumy=0.
            sumt=0.
            do k=1,numdt
              sumt=sumt+iiday(k)
              sumy=sumy+rres(k)
            end do
            avet=sumt/float(numdt)
            avey=sumy/float(numdt)
            sumyt=0.0
            sumtt=0.0
            do k=1,numdt
              sumtt=sumtt+(iiday(k)-avet)**2
              sumyt=sumyt+(rres(k)-avey)*
     &                  (iiday(k)-avet)
            end do
            ratex=sumyt/sumtt
            rate(num+1)=ratex**2
            num=num+1
          end if
        end do
        end do
c  compute MAD of rate (note that rate is in absolute value here
c        chi2(i)=(sumavs/float(num))
        if (num .gt. 2) then
          sumvar=0.0
          do k=1,num
            sumvar=sumvar+rate(k)
          end do
          xrate=sumvar/float(num)
          yrate=xrate
          chi2(i)=xrate
c  remove rates that exceed 2 standard deviations
          rateLimit=Thres*xrate*(1.4826**2)
          kk=0
          do k=1,num
            if (rate(k) .le. rateLimit) then
             kk=kk+1
             rate(kk)=rate(k)
            end if
          end do
            num=kk
          if (num .gt. 2) then 
            sumvar=0.0
            do kk=1,num
              sumvar=sumvar+rate(kk)
            end do

            chi2(i)=sumvar/float(num)
c   note that results of fitting chi2 is dependent on selecting
c     the exponent of num^index ---  1 seems to work best....
            sigavr(i)=chi2(i)/(float(num)**Xscale)
            icount(i)=num
          else
            chi2(i)=1.0e+20
          end if
        else
          chi2(i)=1.0e+20
        end if
      end do
      end if
c  end ((AVRy .eq. 'n') .and. (MADy .eq. 'n')) 
c
c----------------------------------------------------------------------------------------
c
       
c        write(43,*)i,tims(i),chi2(i),num,sigavr(i),yrate,xrate
c      ,(1.4826*xrate)**2

c
c  scan through chi2 values looking for OK values (ie, not 1.0e+20)
c   eliminate bins with bad chi2
c
      kk=0
      do i=1,numav
c        write(44,*)i,tims(i),chi2(i),sigavr(i)
        if (chi2(i) .lt. 0.9e+20) then
          kk=kk+1
          chi2(kk)=chi2(i)
          tims(kk)=tims(i)
          sigavr(kk)=sigavr(i)
        end if
c      write(44,*)i,kk,tims(kk),chi2(kk),sigavr(kk)
      end do
      numav=kk
c
c  Revise sigavr by smoothing.....
c   Fit a power law trend then use that function as sigvar
      do i=1,numav
        sigavr(i)=log(sigavr(i))
        A(i,1)=log(tims(i))
        A(i,2)=1.0
        dum(i)=1.0
        icol(1)=1
        icol(2)=2
      end do
      nmod=2
      call lsqrs(1000,7,numav,nmod,icol,x,e,fit,A,dum,sigavr)
c      print*,"fit to sigvar",(i,x(i),i=1,2)
      do i=1,numav
        sigavr(i)=exp(x(1)*log(tims(i))+x(2))
        if ( OUTy .eq. 'y')
     &    write(52,*)i,tims(i),chi2(i),sigavr(i),icount(i)
      end do
c   Fit a white noise model to AVR
      do i=1,numav
        A(i,1)=1./((tims(i))**3)
        A(i,2)=periodic(182.625,tims(i))
        A(i,3)=periodic(365.25,tims(i))
      end do
      rmmin=1.0d+30
      nmod=3
      do k=1,5
        x(k)=0.0
      end do
      call lsqrs(1000,7,numav,nmod,icol,x,e,fit,A,sigavr,chi2)
      tlen=ent-stt
      sig_scale=365.25*1000.*sqrt(x(1)/(tlen**3))
c      print*,'lsq',x(1),fit,"white noise", 1000*sqrt(x(1)/12.0)," mm",
c     &  "rate error",sig_scale     
      rmmin=fit
      rate_error=sig_scale
      xbest(1)=x(1)
      do i=1,5
        if (i .le. 3) then
          xbest(i)=x(i)
        else
          xbest(i)=0.0
        end if
      end do
      type="white"
      tbest=type
      wnFin=1000.*sqrt(x(1)/(12.))
      do k=2,5
        x(k)=1000.*sqrt(x(k))
      end do
      if ( OUTy .eq. 'y')
     &   write(53,124)type,sig_scale,fit,wnFin,(x(k),k=2,5)
124   format(a6,f10.3,f15.7,2x,5f8.3)

c  seach for optimal power law noise -- two stage grid search
c  first, search from index of 0.5 to 2.5 in 0.1 increments
c  second, starting from optimal index from 1st search, search in 0.01 increments of index
c
      type="pl"
c      print*," "
c      print*,' grid search on power law; coarse search'
      do k=1,4
        icol(k)=k
      end do
      do k=1,21
        plexp=0.5+0.1*float(k-1)
        plexp1=3.0-plexp
        do i=1,numav      
          A(i,4)=1.0/(tims(i)**plexp1)
        end do
        call lsqrsPOS(1000,7,numav,4,icol,x,e,fit,A,sigavr,chi2)
        sig_scale=1000.*365.25*sqrt(x(1)/tlen**3 + x(4)/tlen**plexp1+
     &  x(2)*periodic(182.625,tlen) + x(3)*periodic(365.25,tlen)) 
        wnFin=1000*sqrt(x(1)/12.0)
        fit=fit*sqrt(float(numav-4))/sqrt(float(numav-5))
        if (fit .lt. rmmin ) then
          plexpOpt=plexp
          rmmin=fit
          rate_error=sig_scale
          tbest=type
          do kk=1,5
            xbest(kk)=x(kk)
          end do
        end if
        do kk=2,4
          x(kk)=1000.*sqrt(x(kk))
        end do
        if ( OUTy .eq. 'y')
     &    write(53,124)type,sig_scale,fit,wnFin,(x(kk),kk=2,4),plexp
c
      end do

c  do fine search
c      print*," fine search starting at plexp", plexpOpt
      plexp0=plexpOpt
      do k=-9,9
      if (k .ne. 0) then
        plexp=plexp0 - 0.01*float(k)
        plexp1=3.0-plexp
        do i=1,numav      
          A(i,4)=1.0/(tims(i)**plexp1)
        end do
        call lsqrsPOS(1000,7,numav,4,icol,x,e,fit,A,sigavr,chi2)
        sig_scale=1000.*365.25*sqrt(x(1)/tlen**3 + x(4)/tlen**plexp1 +
     &  x(2)*periodic(182.625,tlen) + x(3)*periodic(365.25,tlen)) 
        wnFin=1000*sqrt(x(1)/12.0)
c        print*,plexp,fit,1000*sqrt(x(1)/12.0)," mm ", sig_scale," mm/yr"
        fit=fit*sqrt(float(numav-4))/sqrt(float(numav-5))
        if (fit .lt. rmmin ) then
          plexpOpt=plexp
          rmmin=fit
          tbest=type
          rate_error=sig_scale
          do kk=1,5
            xbest(kk)=x(kk)
          end do
        end if
        do kk=2,4
          x(kk)=1000.*sqrt(x(kk))
        end do
        if ( OUTy .eq. 'y')
     &    write(53,124)type,sig_scale,fit,wnFin,(x(kk),kk=2,4),plexp
      end if
      end do
c  add flicker and random walk plus seasonal noise; hackl eq 10
      type="flrw"
      do i=1,numav
        A(i,2)=periodic(182.625,tims(i))
        A(i,3)=periodic(365.25,tims(i))
        A(i,4)=1.0/(tims(i)**2)
        A(i,5)=1.0/(tims(i)**1)
        icol(2)=2
        icol(3)=3
        icol(4)=4
        icol(5)=5
c        write(48,*)chi2(i),sigavr(i),(A(i,j),j=1,3)
      end do
      nmod=5
      call lsqrsPOS(1000,7,numav,nmod,icol,x,e,fit,A,sigavr,chi2)
      wnFin=1000*sqrt(x(1)/12.0)
      sig_scale=1000.*365.25*sqrt(x(1)/tlen**3 + x(4)/tlen**2 
     &   + x(5)/tlen + x(2)*periodic(182.625,tlen)
     & + x(3)*periodic(365.25,tlen) )
c note -- to force only a power law estimate of rate, set fit to big number
c    such that flrw results are ignored
      if ( FLRWy .ne. 'y' ) fit=999999.0
c      fit=999999.0
      if (fit .lt. rmmin ) then
        rmmin=fit
        rate_error=sig_scale
        tbest=type
        do k=1,5
          xbest(k)=x(k)
        end do
      end if
      do k=2,5
        x(k)=1000.*sqrt(x(k))
      end do
      if ( OUTy .eq. 'y')
     &  write(53,124)type,sig_scale,fit,wnFin,(x(k),k=2,5)
c
c  compute "fitness"
c
      fit_x=rmmin
       
c  The model prediction
c  
         do j = 1, int(nobs/minav)
            tlen=j*minav
           if ( tbest .eq. "white") then
             pred=xbest(1)*1.0/(tlen**3) 
     &           + xbest(2)*periodic(182.625,tlen)
     &             + xbest(3)*periodic(365.25,tlen)
           end if
           if ( tbest .eq. "flrw" ) then
             pred=xbest(1)*1.0/(tlen**3) 
     &        + xbest(2)*periodic(182.625,tlen)
     &            + xbest(3)*periodic(365.25,tlen)+xbest(4)/(tlen**2)
     &             + xbest(5)/tlen
           end if
           if ( tbest .eq. "pl" ) then
             pred=xbest(1)*1.0/(tlen**3) 
     &        + xbest(2)*periodic(182.625,tlen)
     &            + xbest(3)*periodic(365.25,tlen)
     &            +xbest(4)/(tlen**(3.-plexpOpt))
           end if
         if ( OUTy .eq. 'y')
     &      write(54,123)j,tlen,pred
         end do
123      format(i5,1x,f7.2,e13.4)
      rate_error=rate_error/1000.
      return
      end
      function periodic(T,tau)
c  equation 10 of Hackl
      real*8 tau,periodic
      pi=3.14159265
        per=36.*(T**2)/((pi**2 * tau**4))
        a1=pi*tau/T
        per=per*(sin(a1))**2
        per=per*(sin(a1)/a1 - cos(a1) )**2
        periodic=per
      return
      end
      subroutine lsqrsPOS(md,nd,ic,m,icol,x,e,fit,ar,wt,dr)
c  least squares with crude positivity constraint
c  initially does least squares, then, for x's less than zero
c     removes the columns associated with those x's and redo least squares
c   Subroutine written by John Langbein, Sept 2019
      real*8 d(1113),x(nd), e(nd),ar(md,nd),wt(md)
      real*8 dr(md)
      integer icol(7)
      call lsqrs(md,nd,ic,m,icol,x,e,fit,ar,wt,dr)
c      print*,' lsq ', (i,x(i),i=1,m)," fit",fit

c  check for negative x's -- force to zero and redo lsq
      numx=ic
      ixx=0
      do i=1,m
        if (x(i) .ge. 0 ) then
          ixx=ixx+1
          icol(ixx)=i
        end if
      end do
      
      if (ixx .lt. m) then
        call lsqrs(1000,7,numx,ixx,icol,x,e,fit,ar,wt,dr)
        do k=1,m
          e(k)=0.0
        end do
        do k=1,ixx
          e(icol(k))=x(k)
        end do
        do k=1,m
          x(k)=e(k)
          icol(k)=k
        end do
        fit=fit*sqrt(float(ic-ixx)/float(ic-m))
c        print*,' lsq pos', (i,x(i),i=1,m), " fit", fit
      end if
c  check again for negative x's; this time, punt and set the x <0 to zero
      ixy=0
      do i=1,m
        if (x(i) .lt. 0 ) then
          x(i)=0.0
        else
          ixy=ixy+1
        end if
      end do
      if (ixy .lt. m) then
        sum=0.0
          do j=1,ic
            calc=0.0
            do i=1,m
              calc=calc+ar(j,i)*x(i)
            end do
            sum=sum+((dr(j)-calc)/wt(j))**2
           end do
         fit=sqrt(sum/float(ic-m))
      end if
      return
      end
      subroutine lsqrs(md,nd,ic,m,icol,x,e,fit,ar,wt,dr)
c  least square routine
c mdim row dim of a  max allow num of data
c ndim col dim of a  max allowed mun of model
c ic  num of data
c m num of model para
c icol is an array containing the col # of 'a' for the model
c x  model    e  standard deviation of model
c ws  the data residuals scaled to the aproiri data error
c fit  normalize misfit
c ar  the input 'a' matrix
c wt  the data error
c dr  the data
c covar  the covariance matrix
c     Original subroutine written by John Langbein sometime in the 1980s
      real*8 d(1113),x(nd), e(nd),ar(md,nd),wt(md)
      real*8 dr(md)
      integer icol(7)
      dimension a(1113,157),ev(7),u(1113,7),v(7,7),
     &  ai(7,1113),ws(1113)
      mdim=1113
      ndim=7
c   weight the data and the A matrix
      do 10 i=1,ic
      d(i)=dr(i)/wt(i)
       do 9 j=1,m
9      a(i,j)=ar(i,icol(j))/wt(i)
10    continue
c
c  compute the inverse of A
c
      call svd(mdim,ndim,ic,m,a,ev,u,v,ws,ier)
c examine the eigenvalues
c      write(6,201)(ev(i),i=1,m)
c201   format(//,' the eigenvalues are:',/1x,10(1x,e12.4))
      ie=1
      iflag=0
      do 11 i=2,m
      if (iflag .eq. 1) go to 11
      rat=ev(i)/ev(1)
      if (rat .lt. 5.0e-05) go to 32
      ie=ie+1
      go to 11
32    iflag =1
11    continue
c      write(6,322) ie,m
c322   format(' using',i5,' out of',i5,' eigenvalues')
321   continue
c
      do 15 i=1,m
        do 14 j=1,ic
        at=0.0
          do 13 k=1,ie
13        at=at+v(i,k)*u(j,k)/ev(k)
14      ai(i,j)=at
c      write(7,*)i,(ai(i,k),k=1,ic)
15    continue
c
c compute the model and its error bars
c
      do 17 i=1,m
      x(i)=0.0
      e(i)=0.0
        do 16 j=1,ic
        x(i)=x(i)+ai(i,j)*d(j)
16      e(i)=e(i)+(ai(i,j)**2)
      e(i)=sqrt(e(i))
17    continue

c
c compute the residuals
c
      res=0.0
      do 20 i=1,ic
      ws(i)=0.0
        do 19 j=1,m
19      ws(i)=ws(i)+a(i,j)*x(j)
      ws(i)=d(i)-ws(i)
20    res=res+(ws(i)**2)
      fit=0.0
      if (ic .eq. m) go to 222
      fit=sqrt(res/float(ic-m))
222   continue

      return
      end

c
      subroutine getdata(tx,tz,iday,res,iiday,rres,numdt,nobs,maxm,
     &   nstep,tstep)
c  get data (and time) in interval between tx and tx+dt
c  Subroutine "stolen" from tsfit, originally written for the GAMIT/GLOBK
c    package by Tom Herring
      real*8 res(maxm),rres(maxm)
      integer iday(maxm),iiday(maxm),tstep(100)
      k=0

      do i=1,nobs
        if ((iday(i) .ge. tx) .and. (iday(i) .lt. tz)) then

          do j=1,nstep
            if ((tstep(j) .ge. tx) .and. (tstep(j) .lt. tz)) then
              numdt=0
              return
            end if
          end do

          k=k+1
          iiday(k)=iday(i)
          rres(k)=res(i)
        end if
      end do
      numdt=k
      return
      end
c--------------------------------------------------------------
      subroutine selectpair(m,maxn,tol,t,n,ip,nstep,tstep,twind)
c     Included in MIDAS
c
c     Given a time tag array t(m), select pairs ip(2,n)
c
c     Moves forward in time: for each time tag, pair it with only
c     one future time tag.
c     First attempt to form a pair within tolerance tol of 1 year.
c     If this fails, then find next unused partner.
c     If this fails, cycle through all possible future partners again.
c
c     MIDAS calls this twice -- firstly forward in time, and
c     secondly backward in time with negative tags and data.
c     This ensures a time symmetric solution.

c     2010-10-12: now allow for apriori list of step epochs
c     - do not select pairs that span or include the step epoch

c     input
      integer m, maxn, nstep
      real*8 tol,t(m), tstep(nstep+1)

c     output
      integer n,ip(2,maxn)

c     local
      integer i,j,k,i2,istep
      real*8 dt, fdt,twind
       
      k = 0
      n = 0
      istep = 1
      do i = 1, m
         if( n >= maxn ) exit
         if(t(i) > (t(m)+tol-twind)) exit

c        scroll through steps until next step time is later than epoch 1
         do
            if(istep>nstep) exit
            if(t(i) < tstep(istep)+tol) exit
            istep = istep + 1
         enddo 
         if(istep<=nstep) then
            if(t(i) > (tstep(istep)+tol-twind)) cycle
         endif 

         do j = i+1, m
            if(k<j) k = j
            if(istep<=nstep) then
              if(t(j) > (tstep(istep)-tol)) exit
            endif

            dt = t(j) - t(i)

c         time difference from 1 year (twind years)
            fdt = (dt - twind)

c         keep searching if pair less than one year
            if(fdt<-tol) cycle
c
c         try to find a matching pair within tolerance of 1 year
            if(fdt<tol) then
              i2 = j
            else
c             otherwise, if greater than 1 year, cycle through remaining data
              i2 = k
              dt = t(i2) - t(i)
              if(istep<=nstep) then
                if(t(i2) > (tstep(istep)-tol)) then
                   k = 0
                   cycle
                endif
              endif
              if(k==m) k = 0
              k = k + 1
            endif
c
c         data pair has been found
            n = n + 1
            ip(1,n) = i
            ip(2,n) = i2
            exit
         enddo
      enddo

      end
c
c-----------------------------------------
      subroutine tback(m,t,tb,nstep,tstep,tbstep)
c Included in MIDAS
c     reverses the order of time array t(m)
c     and multiplies result by -1
c     
     
c     input
      integer m, nstep

c     input
      real*8 t(m), tstep(nstep+1), tbstep(nstep+1)

c     output
      real*8 tb(m)

c     local 
      integer i
      
      do i = 1, m
        tb(m+1-i) = -t(i)
      enddo

      do i = 1, nstep
        tbstep(nstep+1-i) = -tstep(i)
      enddo

      end

c---------------------------------------
      real*8 function qmedian(n,arr)
c
c     Wrapper around "qselect" function to find median.
c     This wrapper ensures array "arr" is not overwritten
c     If n is even, qmedian is the mean of the 2 middle numbers
c  Included in MIDAS

      integer n, i, k
      integer, parameter :: maxn=19999 
      real*8 arr(n), work(maxn)
      real*8 qselect
      if (n .le. 1) return
      do i = 1, n
        work(i) = arr(i)
      enddo

      k = (n+1)/2
      qmedian = qselect(k,n,work)

      if (k==n/2) then
         qmedian = 0.5d+0*(qmedian + qselect(k+1,n,work))
      endif

      end

c---------------------------------------
      real*8 function qselect(k,n,arr)
c Included in MIDAS
c
c   Quickselect algorithm by Hoare [1961].
c   Rewritten in modular form to get rid of confusing goto statements.

c   To find the median, set k = (n+1)/2
c   No account is made here as to whether n is odd or even
c
c   Caution: array "arr" gets overwritten (is "in place")
c   This makes if faster for subsequent selects, but care
c   should be taken if "arr" is to be used externally

      integer k, n
      real*8 arr(n)
      integer i, ir, j, l, mid
      real*8 a, temp

      l = 1
      ir = n

      do while (ir > l+1)
         mid = (l+ir)/2
         temp = arr(mid)
         arr(mid) = arr(l+1)
         arr(l+1) = temp

         if (arr(l) > arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
         endif
         if (arr(l+1) > arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
         endif
         if (arr(l) > arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
         endif

         i = l+1
         j = ir
         a = arr(i)

         do
            i = i+1
            if (arr(i) < a) cycle
            do
               j = j-1
               if (arr(j) <= a) exit
            enddo
            if (j < i) exit
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
         enddo

         arr(l+1) = arr(j)
         arr(j) = a
         if (j >= k) ir = j-1
         if (j <= k) l = i
      enddo

      if (ir == 2) then
         if (arr(ir) < arr(l)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
         endif
      endif
      qselect = arr(k)

      end
      subroutine svd(mdim,ndim,m,n,a,w,u,v,temp,ierr)                   
      integer mdim,ndim,m,n,ierr                                        
      real a(mdim,n),w(n),u(mdim,n),v(ndim,n),temp(m)                   
c-----------------------------------------------------------------------
c     this subroutine determines the singular value decomposition       
c          t                                                            
c     a=usv  of a real m-by-n rectangular matrix. it uses the           
c     householder bidiagonalization and the qr algorithms.              
c                                                                       
c     subroutine arguments:                                             
c        mdim must be set to the row dimension of a and u.              
c        ndim must be set to the row dimension of v.                    
c        m is the number of rows of a (and u); m must be >=n.           
c        n is the number of columns of a (and u) and the order of v.    
c        a contains the input matrix to be decomposed; a will not be    
c          altered in svd.                                              
c        w contains the n (non-negative) singular values of a (the      
c           diagonal elements of s) in nonincreasing order. if an error 
c           exit is made, the singular values should be correct         
c           for indices ierr+1,ierr+2,...,n, and they are unordered.    
c        u contains the matrix u (orthogonal column vectors) of         
c           the decomposition. if an error exit is made, the columns    
c           of u corresponding to indices of correct singular values    
c           should be correct.                                          
c        v contains the matrix v (orthogonal) of the decomposition.     
c           if an error exit is made, the columns of v corresponding    
c           to indices of correct singular values should be correct.    
c        temp is a temporary storage array & must be provided by user.  
c        ierr is set to                                                 
c           zero     for normal return,                                 
c           k        if the k-th singular value has not been            
c                    determined after 30 iterations.                    
c-----------------------------------------------------------------------
c     this subroutine is a modified version of the subroutine svd       
c     in the book 'computer methods for mathematical computations,'     
c     by forsythe, malcolm, and moler.                                  
c     modifications by w.h.k. lee & f. luk (1/26/78; 5/30/79).          
c-----------------------------------------------------------------------
      integer i,j,k,l,ii,i1,kk,k1,ll,l1,mn,its                          
      real anorm,f,p,q                                                  
      double precision c,g,h,s,x,y,z,scale                              
      double precision dble,dsqrt,dsign                                 
      ierr = 0                                                          
      do 100 i = 1, m                                                   
         do 100 j = 1, n                                                
  100       u(i,j) = a(i,j)                                             
c--------------- householder reduction to bidiagonal form --------------
      g = 0.0d0                                                         
      scale = 0.0d0                                                     
      anorm = 0.                                                        
      do 300 i = 1, n                                                   
         l = i + 1                                                      
         temp(i) = scale * g                                            
         g = 0.0d0                                                      
         s = 0.0d0                                                      
         scale = 0.0d0                                                  
         if (i .gt. m) go to 210                                        
         do 120 k = i, m                                                
  120       scale = scale + abs(u(k,i))                                 
         if (scale .eq. 0.0d0) go to 210                                
         do 130 k = i, m                                                
            u(k,i) = u(k,i) / scale                                     
  130       s = s + dble(u(k,i))**2                                     
         f = u(i,i)                                                     
         g = -dsign(dsqrt(s),dble(f))                                   
         h = f * g - s                                                  
         u(i,i) = f - g                                                 
         if (i .eq. n) go to 190                                        
         do 150 j = l, n                                                
            s = 0.0d0                                                   
            do 140 k = i, m                                             
  140          s = s + dble(u(k,i)) * dble(u(k,j))                      
            f = s / h                                                   
            do 145 k = i, m                                             
  145          u(k,j) = u(k,j) + f * u(k,i)                             
  150       continue                                                    
  190    do 200 k = i, m                                                
  200       u(k,i) = scale * u(k,i)                                     
  210    w(i) = scale * g                                               
         g = 0.0d0                                                      
         s = 0.0d0                                                      
         scale = 0.0d0                                                  
         if (i .gt. m .or. i .eq. n) go to 290                          
         do 220 k = l, n                                                
  220       scale = scale + abs(u(i,k))                                 
         if (scale .eq. 0.0d0) go to 290                                
         do 230 k = l, n                                                
            u(i,k) = u(i,k) / scale                                     
  230       s = s + dble(u(i,k))**2                                     
         f = u(i,l)                                                     
         g = -dsign(dsqrt(s),dble(f))                                   
         h = f * g - s                                                  
         u(i,l) = f - g                                                 
         do 240 k = l, n                                                
  240       temp(k) = u(i,k) / h                                        
         if (i .eq. m) go to 270                                        
         do 260 j = l, m                                                
            s = 0.0d0                                                   
            do 250 k = l, n                                             
  250          s = s + dble(u(j,k)) * dble(u(i,k))                      
            do 255 k = l, n                                             
  255          u(j,k) = u(j,k) + s * temp(k)                            
  260       continue                                                    
  270    do 280 k = l, n                                                
  280       u(i,k) = scale * u(i,k)                                     
  290    anorm = amax1(anorm,abs(w(i))+abs(temp(i)))                    
  300    continue                                                       
c--------------- accumulation of right-hand transformations ------------
      do 400 ii = 1, n                                                  
         i = n + 1 - ii                                                 
         if (i .eq. n) go to 390                                        
         if (g .eq. 0.0d0) go to 360                                    
         do 320 j = l, n                                                
  320       v(j,i) = (u(i,j) / u(i,l)) / g                              
         do 350 j = l, n                                                
            s = 0.0d0                                                   
            do 340 k = l, n                                             
  340          s = s + dble(u(i,k)) * dble(v(k,j))                      
            do 345 k = l, n                                             
  345          v(k,j) = v(k,j) + s * v(k,i)                             
  350       continue                                                    
  360    do 380 j = l, n                                                
            v(i,j) = 0.                                                 
  380       v(j,i) = 0.                                                 
  390    v(i,i) = 1.                                                    
         g = temp(i)                                                    
         l = i                                                          
  400    continue                                                       
c--------------- accumulation of left-hand transformations -------------
      mn = n                                                            
      if (m .lt. n) mn = m                                              
      do 500 ii = 1, mn                                                 
         i = mn + 1 - ii                                                
         l = i + 1                                                      
         g = w(i)                                                       
         if (i .eq. n) go to 430                                        
         do 420 j = l, n                                                
  420       u(i,j) = 0.                                                 
  430    if (g .eq. 0.0d0) go to 475                                    
         if (i .eq. mn) go to 460                                       
         do 450 j = l, n                                                
            s = 0.0d0                                                   
            do 440 k = l, m                                             
  440       s = s + dble(u(k,i)) * dble(u(k,j))                         
            f = (s / u(i,i)) / g                                        
            do 445 k = i, m                                             
  445          u(k,j) = u(k,j) + f * u(k,i)                             
  450       continue                                                    
  460    do 470 j = i, m                                                
  470       u(j,i) = u(j,i) / g                                         
         go to 490                                                      
  475    do 480 j = i, m                                                
  480       u(j,i) = 0.                                                 
  490    u(i,i) = u(i,i) + 1.                                           
  500    continue                                                       
c--------------- diagonalization of the bidiagonal form ----------------
      do 700 kk = 1, n                                                  
         k1 = n - kk                                                    
         k = k1 + 1                                                     
         its = 0                                                        
c--------------- test for splitting ------------------------------------
  520    do 530 ll = 1, k                                               
            l1 = k - ll                                                 
            l = l1 + 1                                                  
            if ( abs(temp(l))+anorm .eq. anorm) go to 565               
c----- temp(1) is always zero, so there is no exit thru the loop bottom-
  530       continue                                                    
            if ( abs(w(l1))+anorm .eq. anorm) go to 540                 
c--------------- cancellation of temp(l) if l greater than 1 -----------
  540    c = 0.0d0                                                      
         s = 1.0d0                                                      
         do 560 i = l, k                                                
            f = s * temp(i)                                             
            temp(i) = c * temp(i)                                       
            if ( abs(f)+anorm .eq. anorm) go to 565                     
            g = w(i)                                                    
            h = dsqrt( dble(f)*dble(f) + g*g )                          
            w(i) = h                                                    
            c = g / h                                                   
            s = -f / h                                                  
            do 550 j = 1, m                                             
               y = u(j,l1)                                              
               z = u(j,i)                                               
               u(j,l1) = y * c + z * s                                  
  550          u(j,i) = -y * s + z * c                                  
  560       continue                                                    
c--------------- test for convergence ----------------------------------
  565    z = w(k)                                                       
         if (l .eq. k) go to 650                                        
c--------------- shift from bottom 2 by 2 minor ------------------------
c         if (its .eq. 30) go to 1000                                    
c         if (its .eq. 35) go to 1000                                    
         its = its + 1                                                  
         x = w(l)                                                       
         y = w(k1)                                                      
         g = temp(k1)                                                   
         h = temp(k)                                                    
         f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0d0 * h * y)  
         g = dsqrt( dble(f)*dble(f) + 1.0d0 )                           
         f = ((x - z)*(x + z) + h*(y/(f + dsign(g,dble(f))) - h))/x     
c--------------- next qr transformation --------------------------------
         c = 1.0d0                                                      
         s = 1.0d0                                                      
         do 600 i1 = l, k1                                              
            i = i1 + 1                                                  
            g = temp(i)                                                 
            y = w(i)                                                    
            h = s * g                                                   
            g = c * g                                                   
            z = dsqrt( dble(f)*dble(f) + h*h )                          
            temp(i1) = z                                                
            c = f / z                                                   
            s = h / z                                                   
            f = x * c + g * s                                           
            g = -x * s + g * c                                          
            h = y * s                                                   
            y = y * c                                                   
            do 570 j = 1, n                                             
               x = v(j,i1)                                              
               z = v(j,i)                                               
               v(j,i1) = x * c + z * s                                  
  570          v(j,i) = -x * s + z * c                                  
            z = dsqrt( dble(f)*dble(f) + h*h )                          
            w(i1) = z                                                   
c--------------- rotation can be arbitrary if z is zero ----------------
            if (z .eq. 0.0d0) go to 580                                 
            c = f / z                                                   
            s = h / z                                                   
  580       f = c * g + s * y                                           
            x = -s * g + c * y                                          
            do 590 j = 1, m                                             
               y = u(j,i1)                                              
               z = u(j,i)                                               
               u(j,i1) = y * c + z * s                                  
  590          u(j,i) = -y * s + z * c                                  
  600       continue                                                    
         temp(l) = 0.0                                                  
         temp(k) = f                                                    
         w(k) = x                                                       
         go to 520                                                      
c--------------- convergence -------------------------------------------
  650    if (z .ge. 0.0d0) go to 700                                    
c--------------- w(k) is made non-negative -----------------------------
         w(k) = -z                                                      
         do 690 j = 1, n                                                
  690       v(j,k) = -v(j,k)                                            
  700    continue                                                       
c-----arrange singular values in non-increasing order & order u & v also
  805 do 900 k = 1,n                                                    
         p = -1.0                                                       
         do 810 i = k,n                                                 
            if (w(i) .lt. p) go to 810                                  
            p = w(i)                                                    
            j = i                                                       
  810       continue                                                    
         if (j .eq. k) go to 900                                        
         w(j) = w(k)                                                    
         w(k) = p                                                       
         do 820 i = 1,n                                                 
            q = v(i,j)                                                  
            v(i,j) = v(i,k)                                             
  820       v(i,k) = q                                                  
         do 830 i = 1,m                                                 
            q = u(i,j)                                                  
            u(i,j) = u(i,k)                                             
  830       u(i,k) = q                                                  
  900    continue                                                       
      return                                                            
c----- set error -- no convergence to singular value after 30 iterations
 1000 ierr = k                                                          
      return                                                            
      end                                                               

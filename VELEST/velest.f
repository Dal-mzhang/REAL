
c   original velest could be downloaded at http://www.seg.ethz.ch/software/velest.html
c   I met some issues when I complied it on MAC. Here I made some changes to fix them.
c   Miao Zhang (miao.zhang@dal.ca)
c   Main changes by Miao Zhang:
c   1. changed some minor format issue so that it can be complied.
c   2. changed the station 4 char to  6 char
c   3. change the output velocity format (isingle=0),now you can
c      directly use it as input velocity
c   4. fixed a bug (in previous version,the last pick always missing) 
c   
      program velest
c  version 3.3  E.Kissling, Institute of Geophysics, ETH Zurich
c
c                  21. August 1998
c
c  this version is running as: 
c                           
c                               -  velest.f for SUN-Unix 
c                               -  velestUNIXhp.f   for HP-Unix
c
c
c    the program is based on the version 3.0 HRM,UK,EK,and W.Ellsworth
c    of 28.oct93, and on version 3.1 EK 21April95.
c    (see PhD thesis by U.Kradolfer 1989, and by H.R. Maurer 1993
c     at Institute of Gephysics, ETH Zurich, Switzerland)
c
c *********************************************************************
c
c  This version of VELEST can be used to simultaneously invert a large
c  number of earthquake travel time data for hypocentral and velocity
c  model (P and S velocities) parameters. Alternatively (isingle=1) it
c  may also be used as a single event location program.
c
c  This version may read from input data file in the converted (*.cnv)
c  format (ised=0), in a velest-archive type (*arcvel) format (ised=1),
c  and in the new SED (Swiss Seismological Service) format (ised=2).
c  VELEST optionally writes output data files in converted and summary
c  card format and in singleevent mode writes a location file for
c  archive purposes (print output file).
c
c
cEK addition to 3.0:
c  9.nov93  version 3.01 added subr. inputarcvel (ised=1)
c  10.nov93              modified subr. inputparam for model file
c  10.nov93              modified  "        "      for print ouput in
c                                          single event location mode
c  15.nov93              modified in subr. traveltime no calls to 
c                                 routines rejectobs and reviveobs
c  16. nov93             modified in subr. gapcalc to avoid division
c                        zero by zero for call to atan2
c  23. nov93             modified subr inputcnv to recognize end of file
c
c  28. nov93             testing singleevent and multievent mode and
c                        corrections for print output in subr. inputparam
c
c  30. nov93             corrected error in inputparam for isingle and
c                        no of eqs .ne.0:  neqs and legs MUST be one
c                        for isingle=1
c
c implemented dez 1994 EK:   read rotate angle of coordinate system
c                                        as in earlier versions
c                              subr. inputparam read neqs,nshots,rotate
c
c  also dez94 EK:          if ifx(i)=0 normally
c                                   =1 to inhibit adjustement of (rotated)
c                                      y coordinate (see subr. MATRIX)
c                                   (i.e. dtdr(2)=0.0 if ifx(i)=1)
c
c
c EK version 3.3 as of 21August1998:  implemented and improved display of
c                          P and S phase usage and station delays
c                                          
cEK  still to be implemented:  flat Earth approximation
c
cEK
c
c   *vel.out file for STATIS written in subr. STATISLOUT
c
c------------------------------------------------------------------------------
c-----START OF SOURCE-CODE OF PROGRAM VELEST (MAIN PROGRAM) -------------------
c------------------------------------------------------------------------------
c
c
      implicit none
      include 'vel_com.f'
c
      integer i4enabletest
      integer juliam
      integer nittc,istopflag,k,j,i,l,iresflag,ii
      integer ier,iminold,itime
      real damp,xlat,xlon,cpmintot,d1,d2
      real ggti_v(16),gtg_v(16)
      logical better
c
      character cline*80
c
      character*20 ctime
      real cpusec
      cpusec = 0
      ctime = '00:00:00'
c      integer errx,erry,errz
c
c--------------------------------------------
c
c     reset internal cpu-timer:
c
c      call CPUTIMER(cpusec)
c
c test to AVOID runs of velest with compilation-option   NOi4
      i4enabletest=1234567890
c end of this 'mini-test'
c
c     save report of this VELEST-run (write it ev. to file16 in subr. INPUTPARAM
c
c      call DATETIME(ctime)  ! get date&time from system
c      write(headerline(1),'(1x,''>>> Start of program VELEST at '',a20,
c     &              '' <<<'')') ctime
c      write(6,'(1x,''>>> Start of program VELEST at '',a20,
c     &          '' <<<'')') ctime
c      write(6,*)
c
      write(headerline(2),*)'::::::::::::::::::::::::::::::::::::::'
      write(headerline(3),*)' V E L E S T  - Version : January 3, 1995'
c
c      reset...
c
      nitt=0   ! number of iterations made
      nittc=0  !   "          "       for wich CPU is calculated
      nvar=0   ! number of unknowns to solve for
      istopflag=0
c
c     input control-parameters, model and stations
c
      call INPUTPARAM
      if(.not.single_turbo) write(16,*)'~~~ input parameters read.'
c
c     if a file containing all the raypaths (ray-points) of all events
c     is desired, open the file (unit=13):
c
      if(irayout.eq.1) open(13,file=rayfile,
     &                      status='unknown')
      if(idrvout.eq.1) open(21,file=drvfile,
     &                      status='unknown')
c
c     file07 - final hypocenters and travel times ( *.CNV format, for next it.)
c     file11 - summary cards of final hypocenters (for plotting)
c
      if(icnvout.eq.1) open(7,file=cnvfile,
     &                      status='unknown')
      if(ismpout.eq.1) open(11,file=smpfile,
     &                      status='unknown')
      if(isingle.ne.0) open(2,file=velfile,
     &                      status='unknown')  !  statisL-compatible output 
      if(ialeout.eq.1) open(75,file=alefile,
     &                      status='unknown')  !  x,y, ALE of all events
      if(idspout.eq.1) open(76,file=dsprfile,
     &                      status='unknown')  !  x,y, DSPR of all events
      if(irflout.eq.1) open(77,file=rflfile,
     &                      status='unknown')  !  x,y, resid. of reflected ray
      if(irfrout.eq.1) open(78,file=rfrfile,
     &                      status='unknown')  !  x,y, resid. of refracted ray
      if(iresout.eq.1) open(79,file=resfile,
     &                      status='unknown')  !  dist, residual of ray
c
      k=0
ccc      call getsymbol('zirvax')  ! not used in this version
c
 1080 continue     ! START OF THE LOOP OVER  O N E   E V E N T FOR SINGLE-TASK !
c
      if(k.eq.-1)then ! means: end of input-file detected !!!
         close(2)
cek
c    no smp-file for single_event_mode since smp only in SED special format
c    use *.cnv file for plotting etc.
cek
	 if(ismpout.eq.1.and.isingle.eq.1) close(11)
c         call DATETIME(ctime)  ! get date&time from system
c         call CPUTIMER(cpusec)
         if(.not.single_turbo)then
c         write(16,'(x,''>>>   End of program VELEST at '',a20,
c     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
         endif
c         write(6,'(x,''>>>   End of program VELEST at '',a20,
c     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
c         stop'...end...(VELEST was running in SINGLE-EVENT-OPTION)'
      endif
c
c----------------------------------------------------------------------------
c
c     do forward problem first:
c
      nitt=0
c
c     do scaling of damping-factors (and determine ICOUNT):
c
      if(.not.single_turbo)
     &              write(16,*)'~~~ set units for first iteration ...'
c                 1 --> so icount will not be 0 !!!
      call SETUNT(1,invertratio,nsinv,icount,
     &              xythet,stathet,othet,vthet,zthet,scale)
c
c     determine nr of unknowns NVAR & nr of equations LIP to be solved :
c
      if(.not.single_turbo)
     &           write(16,*)'~~~ determine number of unknowns ...'
      call DETNOFUNKNOWNS
c
c     reset G-matrix and RHT-vector :
c
      if(.not.single_turbo) write(16,*)'~~~ reset G and RHT ...'
      do j=1,kvar
         g(j)=0.0
      enddo
      do k=1,nvar
         rht(k)=0.0
      enddo
c
c     for each event do: reset avres&rms, input of the data,
c                        raytracing and accumulate normal equations
c                        ( = setup of matrix G and vector RHT ) :
c
      if(irayout.eq.1) rewind(13) ! rewind file (if desired) with raypaths
      if(idrvout.eq.1) rewind(21) ! rewind file (if desired) with raypaths
      if(.not.single_turbo)
     &        write(16,*)'~~~ input data, raytracing, setup G ...'
      if(.not.single_turbo)then
         write(16,'(/////)')
         write(16,'(''    I N P U T - D A T A  '',/)')
         if(nsp.eq.2) then
         write(16,'(1x,'' eq    origin-time    latitude longitude '',
     &        ''depth    x      y      z    mag ifxOBStot obsP obsS'')')
         else
         write(16,'(1x,'' eq    origin-time    latitude longitude '',
     &                 ''depth    x      y      z    mag ifxOBS'')')
         endif
         endif
      do i=1,legs
         avres(i)=0.
         rms(i)=0.
c        it's the starting: input trial hypocenter and phase-data
         call INPUTDATA(i)
         k=knobs(i)
         if(k.eq.-1)then
            goto 1080   ! means: end of input-file (data) detected
         endif
c
c        for each arrival do raytracing, residual-calc. and setup G & RHT :
c
         do l=1,k
            call TRAVELTIME(i,l,iresflag)
            call SETUPMATRIXG(i,l)
         enddo
         if(isingle.ne.0)then
            if( (knobs(i)-nobswithw0) .lt. nvar )then
               iabort=1
               if(.not.single_turbo)then
                  write(16,*)'knobs(i)-nobswithw0 < nvar !!!'
                  write(16,'('' knobs(i)='',i4,'' nobswithw0='',i4,
     +            '' nvar='',i2)') knobs(i),nobswithw0,nvar
                  write(16,*)'Event cannot be located!!!'
               endif
               write(6,*)'knobs(i)-nobswithw0 < nvar !!!'
               write(6,'('' knobs(i)='',i4,'' nobswithw0='',i4,
     +            '' nvar='',i2)') knobs(i),nobswithw0,nvar
               write(6,*)'Event cannot be located!!!'
               goto 98989
            endif
            if(.not.single_turbo)then
               write(16,*)'~~~ compute singular values of G ...'
            endif
            call SINGULARVALUES(i)
         endif
      enddo
c
c     now all the events are read-in; now check, how many stations
c     appear in the input-data (= nofstainput)
c
      call ACTUALSTATIONS
c
c     save selected columns of g for resolution calculations
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,*)'~~~ store parts of G ...'
      endif
      do k=1,4
         call STOREG(k,k)
      enddo
      if(scale(6).ne.0.0)then
         i=4*neqs+nshot+1
         j=i+nltot-1
         ii=4
         do k=i,j
            ii=ii+1
            call STOREG(k,ii)
         enddo
      endif
c
c     all events: raytracing done, residuals calculated
c                 and normal equations accumulated
c
c     calculate rms & data variance for all events:
c
      if(.not.single_turbo) write(16,*)'~~~ compute RMS and DATVAR ...'
      call RMSDATVAR
c
      if(ittmax.eq.0)then
         if(.not.single_turbo)then
            write(16,*)'WARNING:  ITTmax=0  -->  no iteration is made'
         endif
         write(6,*)'WARNING:  ITTmax=0  -->  no iteration is made'
         goto 9999
      endif
c
c     calculate GAP, average-residuals of ray-types and check,
c     whether variance/msqrd res have decreased (this is
c     always the case if NITT=0 !) or not (in the
c     latter case do hypocenter/model backup):
c
      if(.not.single_turbo) write(16,*)'~~~ save initial datvar ...'
      call CHECKSOLUTION(istopflag,better)
c
c
 10   continue     !  START OF THE LOOP FOR ONE   I T E R A T I O N
c
c     reset CPU-timer for this iteration
c
c      call DATETIME(ctime)  ! get date&time from system
c      call CPUTIMER(cpusec)
      if(isingle.eq.0)then
c         write(6,'(1x,'' finished Iteration #'',i3,''  at'',
c     &             3x,a20,5x,''CPU-sec = '',f7.1)')
c     &             nittc,ctime,cpusec
      endif
      nittc=nittc+1
c
c
      nitt=nitt+1
c
c     print header of this iteration :
c
      cline='----------------------------------------'
     &    //'----------------------------------------'
      if(.not.single_turbo)then
         write(16,'(///,a,///)') cline
         write(16,'(11x,''ITERATION no '',i2,/)') nitt
      endif
c
c     do scaling of damping-factors (and determine ICOUNT) for this iteration:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ set units for this iteration ...'
      endif
      call SETUNT(nitt,invertratio,nsinv,icount,
     &                 xythet,stathet,othet,vthet,zthet,scale)
c
c     set nr of unknowns back to nr for this iteration
c
      if(.not.single_turbo)then
        write(16,*)'~~~ determine number of unknowns for this iter. ...'
      endif
      call DETNOFUNKNOWNS
c
c     print out the damping-factors to be used in this iteration :
c
      if(isingle.eq.0)then
         if(.not.single_turbo)then
            write(16,*)'Damping factors to be used in this iteration:'
            write(16,'(1x,''  Othet='',f8.3,''   XYthet  ='',f8.3,
     &                    ''   Zthet='',f7.3)') othet,xythet,Zthet
            write(16,'(1x,''STAthet='',f8.3,''   Vthet   ='',f8.3,/)')
     &                 stathet,vthet
         endif
      endif
c
c     apply damping to diagonal elements of G-matrix :
c
      if(.not.single_turbo) write(16,*)'~~~ damp G ...'
      call DAMPG
c
c     copy right-hand-side of normal equations (RHT=At*RES) to vector B :
c
      do k=1,nvar
         b(k)=rht(k)
      enddo
c
c
c     Cholesky-decomposition of matrix G :
c
      if(.not.single_turbo) write(16,*)'~~~ solve G * b = RHT ...'
      call LUDECP(g,g,nvar,d1,d2,ier)
      if(ier.ne.0)then   ! matrix not positive definit !!!
         if(.not.single_turbo)then
            write(16,'('' WARNING: error in ludecp     ier='',i5)')ier
         endif
         write(6,'('' WARNING: error in ludecp     ier='',i5)')ier
         STOP'error in (damped) matrix G !!!'
      endif
c
c     solve normal equations (simple elimination; G is now lower triangular) :
c
      call LUELMP(g,b,nvar,b)
c
c     put solution into proper units :
c
      if(.not.single_turbo) write(16,*)'~~~ fix units ...'
      call FIXUNT(b,neqs,nshot,nltot,ksta,scale,
     &            vdamp,itotmodels,inltot,nplay(1))
c
c     check velocity-change;
c     if large, apply step-length-damping to adjustment-vector !
c
      if(.not.single_turbo)then
         write(16,*)'~~~ apply eventually step-length damping ...'
      endif
      call STEPLENGTHDAMP(damp)
c
c     adjust hypocenters, velocities & stationcorrections :
c
      if(.not.single_turbo)then
       write(16,*)'~~~ adjust model (hypocenters, stacorr & veloc.) ...'
      endif
      call ADJUSTMODEL(damp)
c
c     calculate total step-length of this iteration and print it out:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ compute effective step-length ...'
      endif
      call STEPLENGTHCALC
c
c
 1010  continue   ! come here, if backup has been performed   !
c
c     do forward problem for solution of iteration# NITT:
c
c
c     do scaling of damping-factors (and determine ICOUNT) for NEXT iteration:
c
      if(ibackups.lt.4)then
         if(.not.single_turbo)then
            write(16,*)'~~~ set units for next iteration ...'
         endif
         call SETUNT(nitt+1,invertratio,nsinv,icount,
     &                      xythet,stathet,othet,vthet,zthet,scale)
      else
         if(.not.single_turbo)then
            write(16,*)'~~~ set units for this iteration ...'
         endif
         call SETUNT(nitt,invertratio,nsinv,icount,
     &                    xythet,stathet,othet,vthet,zthet,scale)
      endif
c
c     determine nr of unknowns NVAR & nr of equations LIP to be solved for NEXT:
c
      if(.not.single_turbo)then
       if(ibackups.lt.4)then
        write(16,*)'~~~ determine number of unknowns for next iter. ...'
       else
        write(16,*)'~~~ determine number of unknowns for this iter. ...'
       endif
      endif
      call DETNOFUNKNOWNS
c
c     reset G-matrix and RHT-vector :
c
      if(.not.single_turbo) write(16,*)'~~~ reset G and RHT ...'
      do j=1,kvar
         g(j)=0.0
      enddo
      do k=1,nvar
         rht(k)=0.0
      enddo
c
c     for each event do: reset avres&rms, input (if 1st iteration),
c                        raytracing and accumulate normal equations
c                        ( = setup of matrix G and vector RHT ) :
c
      if(irayout.eq.1) rewind(13) ! rewind file (if desired) with raypaths
      if(idrvout.eq.1) rewind(21) ! rewind file (if desired) with raypaths
      if(.not.single_turbo) write(16,*)'~~~ raytracing, setup G ...'
      do i=1,legs
         avres(i)=0.
         rms(i)=0.
         k=knobs(i)
c
c        for each arrival do raytracing, residual-calc. and setup G & RHT :
c
         do l=1,k
            call TRAVELTIME(i,l,iresflag)
            call SETUPMATRIXG(i,l)
c
chrm    Calculate the matrix G for the data resolution matrix
	    if(isingle.ne.0)then
	       gg2(l,1) = 1 * w(l,1)
	       gg2(l,2) = dtdr(1) * w(l,1)
	       gg2(l,3) = dtdr(2) * w(l,1)
	       gg2(l,4) = dtdr(3) * w(l,1)
	    endif
         enddo
c
         if(isingle.ne.0)then
            if( (knobs(i)-nobswithw0) .lt. nvar )then
               iabort=1
               if(.not.single_turbo)then
                  write(16,*)'knobs(i)-nobswithw0 < nvar !!!'
                  write(16,*)'Event cannot be located!!!'
               endif
               write(6,*)'knobs(i)-nobswithw0 < nvar !!!'
               write(6,*)'Event cannot be located!!!'
               goto 98989
            endif
            if(.not.single_turbo)then
               write(16,*)'~~~ compute singular values of G ...'
            endif
            call SINGULARVALUES(i)
         endif
      enddo
c
c     save selected columns of g for resolution calculations
c
c     Store parts of matrix G here (valid for next iteration), before SETUNT
c     for this iteration!
c
      if(.not.single_turbo) write(16,*)'~~~ store parts of G ...'
      do k=1,4
         call STOREG(k,k)
      enddo
      if(scale(6).ne.0.0)then
         i=4*neqs+nshot+1
         j=i+nltot-1
         ii=4
         do k=i,j
            ii=ii+1
            call STOREG(k,ii)
         enddo
      endif
c
c     do scaling of damping-factors (and determine ICOUNT) for CURRENT iter.:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ set units back for current iteration ...'
      endif
      call SETUNT(nitt,invertratio,nsinv,icount,
     &                 xythet,stathet,othet,vthet,zthet,scale)
c
c     set nr of unknowns back to nr for this iteration
c
      if(.not.single_turbo)then
        write(16,*)'~~~ determine number of unknowns for this iter. ...'
      endif
      call DETNOFUNKNOWNS
c
c     all events: raytracing done, residuals calculated
c                 and normal equations accumulated
c
c     calculate rms & data variance for all events:
c
      if(.not.single_turbo)then
         write(16,*)'~~~ compute RMS and DATVAR ...'
      endif
      call RMSDATVAR
c
c     calculate GAP, average-residuals of ray-types and check,
c     whether variance/msqrd res have decreased (this is
c     always the case if NITT=1 !) or not (in the
c     latter case do hypocenter/model backup):
c
      if(.not.single_turbo)then
         write(16,*)'~~~ check solution (better? worse?) ...'
      endif
      call CHECKSOLUTION(istopflag,better)
c
c better=.false. if variance inc'd and residual inc'd, so do readjustments
c if backup was done, but less than 4 times, do forward problem again !!
c
      if(ibackups.lt.4)then
         if(better)then
            if(ibackups.lt.4) ibackups=0
         else
            if(ibackups.eq.0)then
               if(.not.single_turbo)then
                  call NITTOUTPUT(damp)
               endif
            endif
            if(ibackups.lt.4)then
               if(.not.single_turbo)then
                  write(16,'('' msqrd res increased'',
     &                       '' --> hypocenter and model backup.'')')
                  write(16,*)'~~~ solution is worse --> backup ...'
               endif
               call BACKUP
               ibackups=ibackups+1
               call STEPLENGTHCALC
               goto 1010
            endif
         endif
      else
         if(.not.single_turbo)then
            write(16,*)'4 times backup made!'
         endif
      endif
c
c     come here, if no backup made or 4 backups made
c
      if(.not.single_turbo)then
         write(16,*)'~~~ iteration done; output results ...'
      endif
      if(.not.single_turbo) call NITTOUTPUT(damp)
c
      if(.not.single_turbo) write(16,*)'~~~ another iteration? ...'
c
c  IF 4 TIMES 'HALF ADJUSTMENTS MADE', SOLUTION IS FINAL --> STOP ITERATIONS
c
      if(ibackups.eq.4)then
         if(.not.single_turbo)then
            write(16,*)'4 TIMES BACKUP MADE --> sTOP ITERATIONS.'
         endif
         goto 9999
      endif
c
c     IF  STEP-LENGTH < DELMIN ( MINIMUM STEP-LENGTH ) --> STOP ITERATIONS :
c
      if(steplen.gt.0)then  ! otherwise we are during backups...
         if(steplen.lt.delmin)then
            if(.not.single_turbo)then
               write(16,*)'STEP-LENGTH < DELMIN  --> stop iterations.'
            endif
            goto 9999
         endif
      endif
c
c     IF  #ITERATIONS = ITTMAX  --> STOP ITERATIONS :
c
      if(nitt.eq.ittmax)then
         if(.not.single_turbo)then
            write(16,*)'NITT = ITTMAX  --> stop iterations.'
         endif
         goto 9999
      endif
c
c     SINGLE-EVENT-MODE: IF CHANGES IN DATVAR < 1.E-6 --> STOP ITERATIONS :
c
      if(isingle.ne.0.and.istopflag.eq.1)then
         if(.not.single_turbo)then
            write(16,*)'Changes in datvar < 1.e-6  --> stop iterations.'
         endif
         goto 9999
      endif
c
c     DO NEXT ITERATION :
c
      ibackups=0
      GOTO 10
c
c----------------------------------------------------------------------------
c
c     FINAL solution reached!   output final solution of all events,
c                               all residuals and station-corrections:
c
 9999 continue
c
c      call CPUTIMER(cpusec)
      if(isingle.eq.0)then
c        write(6,'(1x,'' finished Iteration #'',i3,''  at'',
c     &             3x,a20,5x,''CPU-sec = '',f7.1)')
c     &             nittc,ctime,cpusec
      endif
      nittc=nittc+1
c
      if(.not.single_turbo)then
         write(16,*)'~~~ final solution reached ...'
         write(16,*)
         if(iresolcalc.gt.0)then
         else
            write(16,*)'~~~ damp G & Chol.-Decomp.: NOT made '//
     &                 '(switch iresolcalc is NOT set !)'
         endif
      endif
      if(iresolcalc.gt.0)then
c
c        G was calculated for next iteration, so compute nvar, kvar etc. for
c        next iteration as well:
         if(.not.single_turbo) 
     &      write(16,*)'~~~ set units for next iteration ...'
         call SETUNT(nitt+1,invertratio,nsinv,icount,
     &                      xythet,stathet,othet,vthet,zthet,scale)
         if(.not.single_turbo) 
     &      write(16,'(a)')' ~~~ determine number of unknowns '//
     &                          'for next iter. ...'
         call DETNOFUNKNOWNS
c
         if(.not.single_turbo) 
     &      write(16,*)'~~~ damp G ...'
         call DAMPG   ! necessary for subr. RESOLCOVAR !!!
         if(.not.single_turbo)then
            write(16,*)'~~~ Cholesky decomposition of G ...'
         endif
         call LUDECP(g,g,nvar,d1,d2,ier)   ! necessary for subr. RESOLCOVAR !
         if(ier.ne.0)then   ! matrix not positive definit !!!
            if(.not.single_turbo)then
              write(16,'('' WARNING: error in ludecp     ier='',i5)')ier
            endif
            write(6,'('' WARNING: error in ludecp     ier='',i5)')ier
            STOP'error in (damped) matrix G !!!'
         endif
      endif
c
      if(isingle.ne.0.and.icoordsystem.eq.2)then
       if(.not.single_turbo)then
          write(16,*)'~~~ compute magnitude ...'
       endif
       do i=1,legs
        if(iyr(i).lt.100)then
           iminold=JULIAM(iyr(i)+1900,imo(i),iday(i),ihr(i),imin(i))
        else
           iminold=JULIAM(iyr(i),imo(i),iday(i),ihr(i),imin(i))
        endif
        call TIMECLEAR(iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &                                                      itime)
        do j=1,knobs(i)
           pt(j,i)=pt(j,i)+(iminold-itime)*60.
        enddo
        call MAGNITUDE(i,itime)
        emag(i)=xmagnitude
c       standard-deviation of emag(i)=sdxmagnitude
c       the magnitudes for each station are now stored on XMAGNI(1...nobs)
       enddo
      endif
c
      if(.not.single_turbo)then
         write(16,*)
         if(isingle.ne.0)then
            write(16,*)'Damping factors to be used in this iteration:'
            write(16,'(1x,''  Othet='',f8.3,''   XYthet  ='',f8.3,
     &                    ''   Zthet='',f8.3)') othet,xythet,Zthet
            write(16,'(1x,''STAthet='',f8.3,''   Vthet   ='',f8.3)')
     &                 stathet,vthet
            write(16,*)
            write(16,*)'Final solution of event#',isingle
            write(16,*)'------------------------------------'
         else !         print final rms-values of all the events
            write(16,59) (rms(i),i=1,legs)
 59         format(' Final rms values: ',10f5.2,/,(19x,10f5.2))
            write(16,'(//)')
         endif
      endif
c
      if(isingle.ne.0)then
         if(icoordsystem.eq.2)then
            call REGION(1,-e(2,1),e(3,1),regionname,nreg,
     &                                   regnamfile,regkoordfile)
         else
            call SDC(e(2,1),e(3,1),xlat,xlon,1) ! calc. LAT/LON
            call REGION(2,xlat,-xlon,regionname,nreg,
     &                       regnamfile,regkoordfile) ! Subr. REGION needs LON E
         endif                                        ! as input (E = positive)
      endif
      if(.not.single_turbo)then
         write(16,*)'~~~ output final hypocenters ...'
      endif
      call FINALHYPOCOUT
      if(.not.single_turbo)then
       if(isingle.ne.0)then
       write(16,'('' Singular values:     '',4f10.4,5x,''ALE ='',f7.3)')
     &                (SV(j),j=1,4) , ale(1)
          write(16,'(1x,''data variance = '',f10.6)') davar1
       endif
      endif
c
      if(.not.single_turbo)then
         write(16,*)'~~~ output final station residuals ...'
         call FINALSTARESI
      endif
      if(iresolcalc.gt.0)then
         if(.not.single_turbo)then
           write(16,*)'~~~ output resolution- & covariance matrices ...'
         endif
         call RESOLCOVAR(davar1)
      endif
c      if (s(2).gt.99.9) then
c         errx = 999
c      else
c      endif
c     if (s(3).gt.99.9) then
c        erry = 999
c      else
c        erry = nint(s(3)*10.)
c      endif
c      if (s(4).gt.99.9) then
c         errz = 999
c      else
c         errz = nint(s(4)*10.)
c      endif
c      write(smpline(46:54),'(3i3)'),errx,erry,errz
cek     if(isingle.eq.0) then
cek         write(11,'(a80)') smpline
cek      endif
      if(ialeout.eq.1)then
         write(75,'(2x,''ALE'',3f10.4,''        .1'')')
     &              -e(2,1),e(3,1), ale(1)
      endif
      if(idspout.eq.1)then
         write(76,'(2x,''DSP'',3f10.4,''        .1'')')
     &              -e(2,1),e(3,1), spread
      endif
c
c     output all the statistics-stuff...
c
      if(isingle.eq.0)then
         write(16,*)'~~~ output statistics ...'
         call STATISTICSOUT
      endif
c
c     close output-files:
c
      if(iturbo.eq.0)then
         close(7)
cek         if(ismpout.eq.1.and.isingle.eq.0) close(11)
      endif
      if(irayout.eq.1)then
         close(13)
         close(21)
      endif
c
c     output cpu-statistics for whole VELEST-run:
c
      if(.not.single_turbo) write(16,*)
c      call CPUTIMER(cpusec)
      if(isingle.eq.0)then
         write(16,*)'CPU - statistics:'
         write(16,*)
      endif
      if(.not.single_turbo)then
         cpmintot=cpusec/60.
         write(16,'(1x,''TOTAL CPU-sec      '',3x,'' ='',f7.1,5x,
     &              ''(CPU-minutes:'',f7.3,'')'')') cpusec,cpmintot
         write(16,*)
      endif
c
      if(isingle.eq.0)then
c         call DATETIME(ctime)  ! get date&time from system
c         write(16,'(x,''>>>   End of program VELEST at '',a20,
c     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
c         write(6,'(1x,''>>>   End of program VELEST at '',a20,
c     &           '' <<<'','' [ CPU-sec : '',f7.1,'']'')') ctime,cpusec
      endif
c
c
c     in case of 'single-event-mode' output NITT, GAP,
c                      SPREAD of resolution and final hypocenter in
c                      statisL-compatible
c
98989 continue   ! return address for 'insufficient data' !!!
      if(isingle.gt.0)then
c         call DATETIME(ctime)  ! get date&time from system
         if(iabort.ne.1)
     &      write(6,'(1x,''GAP = '',i3,''   NITT = '',i2,
     &                ''   D-Spread ='',f5.2)')igap(1),nitt,spread
         write(6,*)
         if(.not.single_turbo)then
            write(16,*)'end of event# ',isingle
            write(16,*)
            write(16,*)
         endif
c
c        output statisL/SED-compatible summary in file02   ( *.VEL )
c
         if(.not.single_turbo)then
         write(16,*)'~~~ output solution in HYPO71-compatible format...'
         endif
c
chrm     Calculation of the data resolution matrix
c        gg2 = Matrix G
c        ggt = G transposed
c	 gtg = G * G transposed
c	 ggti = Inverse of gtg
c        ggg = Generalized inverse
c        drm = data resolution matrix
c
chrm     Setup ggt
c
	 call matrtran(gg2,maxobsperevent,4,ggt)
c
chrm     Setup gtg
c
	 call matrmult(ggt,4,maxobsperevent,gg2,maxobsperevent,4,gtg,4,4)
c
chrm     Damping of gtg ...
c
	 gtg(1,1) = gtg(1,1) + othet
	 gtg(2,2) = gtg(2,2) + xythet
	 gtg(3,3) = gtg(3,3) + xythet
	 gtg(4,4) = gtg(4,4) + zthet
c
chrm     Convert gtg to a vectorized form
c
	 do i=1,4
	    do k=1,4
	       gtg_v((i-1)*4 + k) = gtg(i,k)
	    enddo
	 enddo
	 call matrinv(4,gtg_v,ggti_v)
c
chrm     Reconvert the vectorized matrix ggti_v to a 4 x 4 matrix
c
	 do i=1,4
	    do k=1,4
	       ggti(i,k) = ggti_v((i-1)*4 + k)
	    enddo
	 enddo
c
chrm     Calculate the generalized inverse
c
	 call matrmult(ggti,4,4,ggt,4,maxobsperevent,ggg,4,
     &              maxobsperevent)
c
chrm     Calculate the data resolution matrix
c
	 call matrmult(gg2,maxobsperevent,4,ggg,4,maxobsperevent,drm,
     &             maxobsperevent,maxobsperevent)
c
chrm     Data resolution matrix calculation finished!
c
cek   print output *out.VEL for singleevent locations written by
c     subr. STATISLOUT:
c
         call STATISLOUT
c
c        reset parameters...
c
         iabort=0
         isingle=isingle+1
         nitt=0
         nvar=0
         icount=0
         istopflag=0
         iresflag=0
         nittc=0
         ibackups=0
         ale(1)=0.0
         do i=1,3
            isconstrain(i)=0
         enddo
         iconstrain(1)=0
         emag(1)=0.0
         nmag=0
         do i=1,knobs(1)
            xmagni(i)=0.0
            amx(i)=0.0 ! not necessary in this version;
            prx(i)=0.0 ! just in case, AMX&PRX will be used for *.CNV files also
            istm(i,1)=0
         enddo
         ifixsolution=0
         zmin=zmininput
         goto 1080 ! goto NEXT EVENT (SINGLE-EVENT-MODE !)
      endif
c
c     come here if NOT 'single-event-mode' to terminate program...
c
      close(8) ! close data-input-file ( *.CNV )
c
      goto 99999
c
  998 call OPENERROR('main program','data-input-file FOR008')
c
99999 continue
      END
c
c********       E N D   O F   M A I N   P R O G R A M      ****************
c
cek    begin of vel_io.f
c
      subroutine INPUTPARAM ! old name was: INPUT1 ;
c                           ! this subr. needs another input-format !!
c
c     reads in the control-file (type VELEST.CMN) which contains all the
c     control-parameters for a VELEST-run
c
      implicit none
      include 'vel_com.f'
      real olat,olon,avelo,z,dx,dy,dz,rotate
      integer ifil,j,ifl,mode,icc,ml,jndex,k
c
      integer trimlen, i,n
      character*1 reflch, cns,cew
      character*40 titl
      character*80 card,line(32), titleline
      logical lexist
c
c     Open control-inputfile:
c
      inquire(file='velest.cmn',exist=lexist)
      if(.not.lexist)then
         stop'INPUTPARAM>>> control file `velest.cmn` not found!'
      endif
cVMS      open(10,file='velest.cmn',status='old',err=9910,readonly)
      open(10,file='velest.cmn',status='unknown',err=9910)
c
c     input center of coordinate system
c
      i=0
 111  read(10,'(a)',end=222) card
      if(card(1:1).eq.'*') goto 111
      i=i+1
      if(i.gt.32) stop'INPUTPARAM>>> control-file not correct!'
      line(i)=card
      goto 111
 222  continue
      if(.NOT.(i.eq.24.or.i.eq.32.))
     &    stop'INPUTPARAM>>> control-file not correct!'
      read(line(1),'(a)') titleline
      read(line(2),*) olat,olon,icoordsystem,zshift,itrial,ztrial,ised
      read(line(3),*) neqs,nshot,rotate
      read(line(4),*) isingle, iresolcalc
      read(line(5),*) dmax,itopo,zmininput,veladj,zadj,lowveloclay
      read(line(6),*) nsp, swtfac,vpvs, nmod
      read(line(7),*) othet,xythet,zthet,vthet,stathet
      read(line(8),*) nsinv,nshcor,nshfix, iuseelev,iusestacorr
      read(line(9),*) iturbo, icnvout,istaout,ismpout
      read(line(10),*) irayout,idrvout,ialeout,idspout,
     &                irflout,irfrout,iresout
      read(line(11),*) delmin,ittmax,invertratio
      read(line(12),'(a)') modelfilename
      read(line(13),'(a)') stationfilename
      read(line(14),'(a)') seismofilename
      read(line(15),'(a)') regnamfile
      read(line(16),'(a)') regkoordfile
      read(line(17),'(a)') topo1file
      read(line(18),'(a)') topo2file
      read(line(19),'(a)') phasefile
      read(line(20),'(a)') shotfile
c     Output files:
      read(line(21),'(a)') outfile
      if(outfile.eq.' ') outfile='vel.out'
      read(line(22),'(a)') velfile
      if(velfile.eq.' ') velfile='velout.vel'
      read(line(23),'(a)') cnvfile
      if(cnvfile.eq.' ') cnvfile='velout.cnv'
      read(line(24),'(a)') stafile
      if(stafile.eq.' ') stafile='velout.sta'
c
c     the next few ouputfiles are not very often used, therefore
c     either all or none of them have to be specified in the controlfile:
c
      if(i.eq.32)then
         read(line(25),'(a)') smpfile
         if(smpfile.eq.' ') smpfile='velout.smp'
         read(line(26),'(a)') rayfile
         if(rayfile.eq.' ') rayfile='velout.ray'
         read(line(27),'(a)') drvfile
         if(drvfile.eq.' ') drvfile='velout.drv'
         read(line(28),'(a)') alefile
         if(alefile.eq.' ') alefile='velout.ale'
         read(line(29),'(a)') dsprfile
         if(dsprfile.eq.' ') dsprfile='velout.dspr'
         read(line(30),'(a)') rflfile
         if(rflfile.eq.' ') rflfile='velout.rfl'
         read(line(31),'(a)') rfrfile
         if(rfrfile.eq.' ') rfrfile='velout.rfr'
         read(line(32),'(a)') resfile
         if(resfile.eq.' ') resfile='velout.res'
      endif
c
      single_turbo=.false.
      if(isingle.eq.1.and.iturbo.eq.1) single_turbo=.true.
c
c     open the main-output-file:
c
      if(.not.single_turbo)then
cVMS         open(16,file=outfile,status='new',carriagecontrol='list')
         open(16,file=outfile,status='unknown')
         write(16,'(a)') headerline(1)
         write(16,*)
         write(16,'(a)') headerline(2)
         write(16,'(a)') headerline(3)
         write(headerline(3),*)' (Authors: see source code)'
         write(16,'(a)') headerline(3)
         write(16,'(a)') headerline(2)
         write(16,*)
         write(16,*)
         write(16,*)'Title of this VELEST run:'
         write(16,*)
         n=trimlen(titleline)
         write(16,*) titleline(1:n)
         do i=1,n
            titleline(i:i)='-'
         enddo
         write(16,*) titleline(1:n)
         write(16,*)
         write(16,*)
         write(16,*)'Current array-dimensions of program VELEST:'
         write(16,*)
         write(16,'('' Max. number of '',
     &              ''EARTHQUAKES for simult. inversion IEQ = '',i3)')
     &              ieq
         write(16,'('' Max. number of '',
     &              ''SHOTS for simult. inversion    INSHOT = '',i3)')
     &              inshot
         write(16,'('' Max. number of '',
     &              ''(different) 1D-MODELS      ITOTMODELS = '',i3)')
     &              itotmodels
         write(16,'('' Max. number of '',
     &              ''LAYERS per one 1D-model        INLTOT = '',i3)')
     &              inltot
         write(16,'('' Max. number of '',
     &              ''STATIONS in stationlist           IST = '',i3)')
     &              ist
         write(16,'('' Max. number of '',
     &              ''OBSERVATIONS per event MAXOBSPEREVENT = '',i3)')
     &              maxobsperevent
         write(16,*)
         write(16,*)
      endif
c
cek 210495: check for reasonable switch combinations nmod,nsp
c
      if(nsp.lt.1) nsp=1
      if(nsp.gt.3) nsp=3
c nsp=1 P-model and p-data only
c nsp=2 P-model (+data) and S-model (+data) independently used
c nsp=3 P-model (+data) and S-model dependent on P-model (vpvs-factor used)
      if(nsp.eq.2) then
        nmod=2
      else
        nmod=1
      endif
c
cek 210495
c
c     print the input-parameters
c
      if(.not.single_turbo)then
         write(16,'(///)')
         write(16,*)'   INPUT - P A R A M E T E R S :'
         write(16,'(//)')
      write(16,'(a)')'***'
      write(16,'(a)')'***  olat       olon   icoordsystem      '//
     &              'zshift   itrial ztrial    ised'
      write(16,'(1x,f10.6,f11.6,6x,i1,10x,f7.3,7x,i1,3x,f6.2,6x,i1)')
     &          olat,olon,icoordsystem,zshift,itrial,ztrial,ised
      write(16,'(a)')'***'
      write(16,'(a)')'*** neqs   nshot    rotate'
      write(16,'(1x,i6,1x,i6,3x,f6.1)') neqs,nshot,rotate
      write(16,'(a)')'***'
      write(16,'(a)')'*** isingle   iresolcalc'
      write(16,'(7x,i1,10x,i1)') isingle, iresolcalc
      write(16,'(a)')'***'
      write(16,'(a)')'*** dmax    itopo    zmin'//
     &              '     veladj    zadj   lowveloclay'
      write(16,'(2x,f7.2,5x,i1,4x,f6.2,6x,f5.2,3x,f5.2,7x,i1)')
     &          dmax,itopo,zmininput,veladj,zadj,lowveloclay
      write(16,'(a)')'***'
      write(16,'(a)')'*** nsp    swtfac   vpvs       nmod'
      write(16,'(5x,i1,5x,f5.2,3x,f6.3,7x,i2)') nsp, swtfac, vpvs, nmod
      write(16,'(a)')'***'
      write(16,'(a)')'***   othet   xythet    zthet    vthet   stathet'
      write(16,'(4x,5(f7.2,2x))') othet,xythet,zthet,vthet,stathet
      write(16,'(a)')'***'
      write(16,'(a)')'*** nsinv   nshcor   nshfix'//
     &              '     iuseelev    iusestacorr'
      write(16,'(7x,3(i1,7x),4x,i1,12x,i1)') nsinv, nshcor, nshfix,
     &                         iuseelev,iusestacorr
      write(16,'(a)')'***'
      write(16,'(a)')'*** iturbo    icnvout   istaout   ismpout'
      write(16,'(7x,i1,9x,i1,9x,i1,9x,i1)') iturbo,
     &                                      icnvout,istaout,ismpout
      write(16,'(a)')'***'
      write(16,'(a)')'*** irayout   idrvout   ialeout   idspout'//
     &               '   irflout   irfrout   iresout'
      write(16,'(7x,i1,6(9x,i1))') irayout,idrvout,ialeout,idspout,
     &                            irflout,irfrout,iresout
      write(16,'(a)')'***'
      write(16,'(a)')'*** delmin   ittmax   invertratio'
      write(16,'(3x,f6.3,6x,i2,9x,i2)') delmin, ittmax, invertratio
      write(16,'(a)')'***'
      write(16,'(a)')'*** Modelfile:'
      write(16,'(a)') modelfilename(1:trimlen(modelfilename))
      write(16,'(a)')'***'
      write(16,'(a)')'*** Stationfile:'
      write(16,'(a)') stationfilename(1:trimlen(stationfilename))
      write(16,'(a)')'***'
      write(16,'(a)')'*** Seismofile:'
      write(16,'(a)') seismofilename(1:trimlen(seismofilename))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with region names:'
      write(16,'(a)') regnamfile(1:trimlen(regnamfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with region coordinates:'
      write(16,'(a)') regkoordfile(1:trimlen(regkoordfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File #1 with topo data:'
      write(16,'(a)') topo1file(1:trimlen(topo1file))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File #2 with topo data:'
      write(16,'(a)') topo2file(1:trimlen(topo2file))
      write(16,'(a)')'***'
      write(16,'(a)')'*** DATA INPUT files:'
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with Earthquake data:'
      write(16,'(a)') phasefile(1:trimlen(phasefile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with Shot data:'
      write(16,'(a)') shotfile(1:trimlen(shotfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** OUTPUT files:'
      write(16,'(a)')'***'
      write(16,'(a)')'*** Main print output file:'
      write(16,'(a)') outfile(1:trimlen(outfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with single event locations:'
      write(16,'(a)') velfile(1:trimlen(velfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with final hypocenters in *.cnv format:'
      write(16,'(a)') cnvfile(1:trimlen(cnvfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with new station corrections:'
      write(16,'(a)') stafile(1:trimlen(stafile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with summary cards (e.g. for plotting):'
      write(16,'(a)') smpfile(1:trimlen(smpfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with raypoints:'
      write(16,'(a)') rayfile(1:trimlen(rayfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with derivatives:'
      write(16,'(a)') drvfile(1:trimlen(drvfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with ALEs:'
      write(16,'(a)') alefile(1:trimlen(alefile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with Dirichlet spreads:'
      write(16,'(a)') dsprfile(1:trimlen(dsprfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with reflection points:'
      write(16,'(a)') rflfile(1:trimlen(rflfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with refraction points:'
      write(16,'(a)') rfrfile(1:trimlen(rfrfile))
      write(16,'(a)')'***'
      write(16,'(a)')'*** File with residuals:'
      write(16,'(a)') resfile(1:trimlen(resfile))
      write(16,'(a)')'***'
      write(16,'(///)')
      endif
c
c     remarks concerning OLAT & OLON:
c     Both olat and olon are input in decimal-degrees.
c     For LAT N(orth) and LON W(est) their values are positive.
c     For LAT S(outh) and LON E(ast) their values must be negative!
c
c     For more details see comments on top of this program-source!
c
      zmin=zmininput
c
      if(.not.single_turbo)then
         write(16,*)'Origin of cartesian coordinates:'
         if(olon.le.0.0)then
           write(16,'(6x,f10.5,2hN ,f11.5,1hE)')olat,-olon
         else
           write(16,'(6x,f10.5,2hN ,f11.5,1hW)')olat,olon
         endif
         if(olon.eq.0.0.and.olat.eq.0.0)then
            write(16,*)'Origin BERNE is taken!'
            write(16,*)'OLAT=46.95240 n OLON=7.439583 e'
         endif
         write(16,*)
         write(16,*)' X,Y-AXES rotated clockwise from North'
         write(16,*)'          Rotation angle (rotate):'
         write(16,'(12x,7hrotate=,1x,f6.1,1x,5hdegr.)')rotate
         write(16,*)
c
         if(itrial.gt.0)then
            write(16,*)
            write(16,*)' ============================================'
            write(16,*)' trial epicenter ~ earliest station'
            write(16,*)' trial depth = ztrial = ',ztrial
            write(16,*)' ============================================'
            write(16,*)
         endif
c
         if(icoordsystem.eq.2)then
            write(16,*)'SWISS COORDINATES will be used instead of the',
     &                 ' short distance conversion!'
            write(16,*)'origin BERNE: x=600.km, y=200.km'
            write(16,*)'x is positive towards east,',
     &                 ' y is positive towards north'
ccc            write(16,*)'---> no rotation can be performed !!!'
         else
            write(16,*)'icoordsystem = ',icoordsystem
            write(16,*)'normal SHORT DISTANCE CONVERSION will be made'
            write(16,*)'x is positive towards west,',
     &                 ' y is positive towards north'
         endif
      endif
c
c     use setorg to set up coordinate system
c
      ifil=0
      if(.not.single_turbo) ifil=16
c
      if(icoordsystem.ne.2) call SETORG(olat,olon,rotate,ifil)
c
c  now origin is setup.
c
4     format(f5.2,i5,4f5.2,f6.2,3i1,5f6.2)
c
      if(dmax.eq.0.)then
         if(.not.single_turbo)then
            write(16,*)'WARNING: dmax was zero ! ... set to 150 km !'
         endif
         dmax=150.
      endif
c
      rmsmin=0.0
      if(invertratio.le.0) invertratio=999  ! do not invert for sta-corr & model
      if(.not.single_turbo)then
         if(iuseelev.eq.0)then
            write(16,*)
            write(16,*)'Station-elevations internally set to ZERO !'
            write(16,*)'(but correctly printed in file12)'
            write(16,*)
         endif
         if(iusestacorr.eq.0)then
            write(16,*)
            write(16,*)'Station-corrections set to ZERO !'
            write(16,*)'(if you do NOT invert for station-corrections,'
            write(16,*)' then these 0.0-values are printed in file12 !)'
            write(16,*)
         endif
      endif
c***  nsinv=  0 no inversion for station corrections
c***          1 inversion for station corrections
c
c***  legs is the total number of shots and quakes
cek
cek   legs and neqs MUST be 1 in single_event mode
c
      if(isingle.ne.1) then
         legs=neqs+nshot
      else
         legs=1
         neqs=1
         nshot=0
      endif
c
      if(.not.single_turbo)then
         write(16,'(///)')
         write(16,*)'   INPUT - M O D E L :'
         write(16,'(//)')
         write(16,*)'Model(s) read from file :'
         write(16,*) modelfilename
         write(16,*)
      endif
      close(10)
c
c   ***********
c NOW READ model file
c   ***********
c
cVMS      open(10,file=modelfilename,status='old',err=9911,readonly)
      open(10,file=modelfilename,status='old',err=9911)
c
      ireflector=0
cek
      if(nsp.eq.3) then
         write(16,7788)
7788  format(//,2x,'attention nsp=3! >s-data used for p-model too!',//)
      endif
c
cek  modifications for new model file
      read(10,'(a40)') titl
      if(.not.single_turbo) write(16,'(a)') 
     +                 ' model file title: ',titl
cek      read(10,*) (nplay(j),j=1,nmod)
8     format(10i5)
      do 14 i=1,nmod
      if(.not.single_turbo)then
         write(16,17) i
17       format(1x,'Velocity structure for model',i4,' :')
         write(16,*)
         write(16,11)
11       format(1h ,'layer    vel   depth   vdamp  reflector')
      endif
      titl=' '
      read(10,'(i3)') nplay(i)
      do 9 j=1,nplay(i)
       if(j.eq.1)then
         read(10,1212) vp(i,j),hp(i,j),vdamp(i,j),reflch,titl
1212     format(f5.2,5x,f7.2,2x,f7.3,3x,a1,1x,a40)
       else
         read(10,12) vp(i,j),hp(i,j),vdamp(i,j),reflch
12       format(f5.2,5x,f7.2,2x,f7.3,3x,a1)
       endif
      if(reflch.ne.' ')then
         if(reflch.eq.'m'.or.reflch.eq.'M')then
            reflchar=reflch
            ireflector=j
            if(.not.single_turbo)then
               if(vp(i,ireflector).gt.8.0)then
                  write(16,*)'WARNING ::: velocity ABOVE reflector is',
     &                       ' greater than 8.0 km/s  .....'
                  write(6,*)'WARNING ::: velocity ABOVE reflector is',
     &                       ' greater than 8.0 km/s  .....'
                  write(6,*)
               endif
            endif
         else
            if(.not.single_turbo)then
              write(16,*)'WARNING:'
              write(16,*)'Reflector indicated in velocity-model is',
     &                   ' marked with a : ',reflch
              write(16,*)'Only   m   or   M   are allowed!'
              write(16,*)'Reflected phases will be ignored (wrong mark)'
c              write(6,*)'Reflected phases will be ignored (wrong mark)'
c              write(6,*)
            endif
         endif
      endif
      if(j.eq.1)then
         if(.not.single_turbo)then
            write(16,'(2x,i2,3x,2f7.2,1x,f7.2,1x,a1,1x,a40)') j,
     &                vp(i,j),hp(i,j),vdamp(i,j),reflch,titl
         endif
         titl=' '
      else
         if(.not.single_turbo)then
            write(16,13) j,vp(i,j),hp(i,j),vdamp(i,j),reflch
13          format(2x,i2,3x,2f7.2,1x,f7.2,1x,a1)
         endif
      endif
c
      if(j.gt.1.and.hp(i,j).lt.0.0)then
         write(6,*)'WARNING:'
         write(6,*)'Only top of the first layer can be negative !!'
         write(6,*)
         if(.not.single_turbo)then
            write(16,*)'WARNING:'
            write(16,*)'Only top of the first layer can be negative !!'
         endif
      endif
c
    9 continue
c
c    calculate and print average velocities of the model i :
c
      if(.not.single_turbo)then
        ifl=1
        write(16,*)
        write(16,*)'Calculation of average velocity starts at layer # ',
     &             ifl
        avelo=0
        do j=ifl+1,nplay(i)
           avelo=avelo + ( hp(i,j)-hp(i,j-1) ) * vp(i,j-1)
           write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &               ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,j-1),hp(i,j),vp(i,j-1),avelo/(hp(i,j)-hp(i,ifl)),
     &           hp(i,j)
        enddo
        write(16,*)
c
        ifl=2
        write(16,*)
        write(16,*)'Calculation of average velocity starts at layer # ',
     &             ifl
        avelo=0
        do j=ifl+1,nplay(i)
           avelo=avelo + ( hp(i,j)-hp(i,j-1) ) * vp(i,j-1)
           write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &               ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,j-1),hp(i,j),vp(i,j-1),avelo/(hp(i,j)-hp(i,ifl)),
     &           hp(i,j)
        enddo
        write(16,*)
        write(16,*)
      endif
c
14    continue
      i f (ireflector.ne.0) t h e n
         if(.not.single_turbo)then
         write(16,*)
        write(16,'(1x,''Phases in the input-phaselist marked with a  '',
     &            a1,/,'' are treated as reflections from the bottom '',
     &            ''of layer nr. '',i2)') reflchar,ireflector
            if(lowveloclay.eq.1)then
               write(16,*)'Switch LOWVELOCLAY is set to 1, but'//
     &                    ' reflected phases are allowed to occur'
               write(16,*)'This is improper (no low velocity-layers '//
     &                    'allowed for reflected waves!!)'
               lowveloclay=0
               write(16,*)'LOWVELOCLAY now set to 0 '
               write(6,*)'Switch LOWVELOCLAY is set to 1, but'//
     &                    ' reflected phases are allowed to occur'
               write(6,*)'This is improper (no low velocity-layers '//
     &                    'allowed for reflected waves!!)'
               write(6,*)'LOWVELOCLAY now set to 0 '
            endif
         endif
      e n d i f
c
c
c   ***********
c NOW READ station file
c   ***********
c
c    read in station data:
c
      if(.not.single_turbo)then
         write(16,'(///)')
         write(16,*)'   INPUT - S T A T I O N S :'
         write(16,'(//)')
         write(16,*)'Station-parameters read from file :'
         write(16,*) stationfilename
         write(16,*)
      endif
      close(10)
cVMS      open(10,file=stationfilename,status='old',err=9912,readonly)
      open(10,file=stationfilename,status='old',err=9912)
      read(10,1) fm
1     format(a80)
      if(.not.single_turbo)then
         write(16,*)
         write(16,2)
2        format(1h ,3x,' stn latitude longitude elev',
     &   3x,'x      y      z',5x,'ptcor stcor model icc')
      endif
      nsta=0
c
10    nsta=nsta+1
      read(10,fm) stn(nsta),xla(nsta),cns,xlo(nsta),cew,
     &            ielev(nsta),mode,icc,
     &            ptcor(nsta),stcor(nsta)
      call CASEFOLD(cns)
      call CASEFOLD(cew)
      if(cns.eq.'S') xla(nsta)=-xla(nsta)
      if(cew.eq.'E') xlo(nsta)=-xlo(nsta)
c
      if(stn(nsta).eq.' ') goto 41
cek      if(mode.eq.0) mode=1
      mode=1
c
      if(nsp.eq.2) then
         model(2*nsta-1)=mode
         model(2*nsta)=mode+1
      else
         model(nsta)=mode
      endif	 
      z=-ielev(nsta)/1000.
      if(icoordsystem.eq.2)then
         call GEOKO(dx,dy,xla(nsta),-xlo(nsta),-1) ! calc. cart. coord.
         dx=-dx
      else
         call SDC(dx,dy,xla(nsta),xlo(nsta),-1) ! calc. cart. coord.
      endif
      dz=z
      x(nsta,1)=dx
      x(nsta,2)=dy
      if(iuseelev.eq.1)then
         x(nsta,3)=dz
      else
         x(nsta,3)=0.     ! station-elevation not used !
      endif
      if(iusestacorr.eq.0)then
         ptcor(nsta)=0.0  ! initial (from input) station-corrections not used !
         stcor(nsta)=0.0  !   "          "               "            "    "
      endif
      if(.not.single_turbo)then
         if(icoordsystem.eq.2) dx=-dx
         if(cns.eq.'S') xla(nsta)=-xla(nsta)
         if(cew.eq.'E') xlo(nsta)=-xlo(nsta)
         write(16,27) nsta,stn(nsta),xla(nsta),cns,xlo(nsta),cew,
     &                ielev(nsta),dx,dy,dz,
     &                ptcor(nsta),stcor(nsta),mode,icc
27       format(1x,i3,1x,a6,f7.4,a1,f8.4,a1,1x,i5,3f7.2,1x,2f6.2,2i4)
         if(cns.eq.'S') xla(nsta)=-xla(nsta)
         if(cew.eq.'E') xlo(nsta)=-xlo(nsta)
      endif
c
c   an icc of 0 holds the station delay fixed.
c   the highest nonzero icc station has its p delay
c   held fixed but its s delay is allowed to float.
c
      map1(nsta)=icc
      goto 10
41    continue
      close(10)  ! stationfile
c
c     read in seismo data (seismometer-specifications and so on):
c
      if(.not.single_turbo)then
         if(seismofilename.ne.' ')then
            write(16,*)
            write(16,*)'Seismo-parameters read from file :'
            write(16,*) seismofilename
            write(16,*)
            write(16,*)'This file contains station-informations which'
            write(16,*)'are used for magnitude-determination'
            write(16,*)'(seismometer-constant,filterparameters, etc.)'
            write(16,*)
            write(16,*)
         else
            write(16,*)
            write(16,*)'NO Seismo-parameter-file specified'
            write(16,*)'--> NO magnitudes calculated!'
            write(16,*)
            write(16,*)
            write(16,*)
         endif
      endif
      if(seismofilename.ne.' ')
cVMS     &open(10,file=seismofilename,status='old',err=9913,readonly)
     &open(10,file=seismofilename,status='old',err=9913)
c     the seismo-parameters will be read in later; we must know the eventtime
c     before we can pick the time-dependent (!) seismo-parameters !!!
c
c***  nsta is the number of stations read in
      nsta=nsta-1
c
      call MAXII(nsta,map1,ml,jndex) ! determine MAX (=ml) of map1 (icc-values)
c***  ksta is the number of total station corrections
c***    to invert for
      ksta=0
      if(nsinv.ne.0)then  !  invert for station corrections
c                            the last station for p  will be held fixed.
         ksta=ml-1
         if(nsp.eq.2) ksta=2*ksta + 1
      endif
c
c    find the total number of layers
c
c**  laysum(i) is an index for the location of the first
c***  layer of model i
      laysum(1)=1
      laysum(2)=nplay(1)+1
      if(nmod.lt.3) goto 50
      do 42 i=3,nmod
42    laysum(i)=laysum(i-1)+nplay(i-1)
50    continue
c***  nltot is the total number of velocity layers
c***  to invert for
      nltot=0
      do 18 i=1,nmod
      nltot=nltot+nplay(i)
18    continue
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,19) neqs,nshot,nltot
      endif
 1619 format(10x,3i20)
19    format(3x,'neqs=',i5,'  nshot=',i5,'  nltot=',i5,/)
c
      do 25 i=1,itotmodels     ! max. number of 1D-models
      do 25 j=1,inltot         ! max. number of layers per model
25    thkp(i,j)=0.0
      do 22 i=1,nmod
      if(nplay(i).eq.1) goto 22
      k=nplay(i)-1
      do 23 j=1,k
23    thkp(i,j)=hp(i,j+1)-hp(i,j)
22    continue
      RETURN
c
 9910 call OPENERROR('inputparam','control-input-file FOR010')
      return
 9911 call OPENERROR('inputparam','model-input-file (FOR010)')
      return
 9912 call OPENERROR('inputparam','station-input-file (FOR010)')
      return
 9913 call OPENERROR('inputparam','seismo-input-file (FOR010)')
      return
c
      end ! of subr. inputparam
c
c
      subroutine INPUTDATA(i)  ! old name:  INPUT2
c
c     reads in the phaselists with the earthquake and/or shot data
c
      implicit none
      include 'vel_com.f'
c
      integer     i
      integer     nobs,ntobs,jj,nobsread,i1,i2,i3,icc,jjmin1,jk,
     &            jjmin,jshot,j,k,iunit,ll,itest,ie,iphaseteststopflag,
     &            l,nobsp,nobss
      real  wsum,xlat,xlon,alon,depth,z,ss1
      real  xxlat,aalon,xxx
      character*1 sc,ss,cphase(maxobsperevent)
      character*1 rmk1(maxobsperevent), rmk2(maxobsperevent),
     &            eventtype, cns,cew
      character*6 sta(maxobsperevent)
      integer ipwt(maxobsperevent)
      real sec(maxobsperevent)
c
      SAVE iphaseteststopflag  ! forget never that error was detected!
c
      data sc,ss/'s','S'/
c
c     input phase list and set-up initial trial hypocenter
      iunit=8
      if(i.eq.1)then
         if(neqs.gt.0) open(iunit,file=phasefile,
cVMS     &                      status='old',err=998,readonly)
     &                      status='old',err=998)
         if(nshot.gt.0) open(9,file=shotfile,
cVMS     &                       status='old',err=999,readonly)
     &                       status='old',err=999)
      endif
cek
cek      re-entry for reading next event if first event has less
cek      than three observations
10001 continue
cek
      nobs=0
      nobswithw0=0
      ntobs=0
      wsum= 0.0
c
      if(i.le.neqs)then
         iunit=8    ! EQS
      else
         iunit=9    ! SHOTS
      endif
c
      xlat=0.0
      xlon=0.0
      do jj=1,maxobsperevent
         cphase(jj)=' '
      enddo
c
cek  i1-switch is dummy
      i1=0
c
cek     EQS:   i1=ifx(i)     ! no longer in use !!!
cek     SHOTS: i1=icc
c
      if(ised.eq.0)then
         call INPUTCNV(iunit,nobsread,
     &                 sta,iyr(i),imo(i),iday(i),ihr(i),imin(i),sec,
     &                 rmk1,rmk2,cphase,ipwt,amx,prx,
     &                 xlat,alon,emag(i),depth,e(1,i),
     &                 i1,i2,i3,eventtype)
      endif
      if(ised.eq.1)then
         call INPUTARCVEL(iunit,nobsread,
     &                 sta,iyr(i),imo(i),iday(i),ihr(i),imin(i),sec,
     &                 rmk1,rmk2,cphase,ipwt,amx,prx,
     &                 xlat,alon,emag(i),depth,e(1,i),
     &                 i1,i2,i3,eventtype)
      endif
      if(ised.eq.2)then
         call INPUTSED_NEW(iunit,nobsread,
     &                 sta,iyr(i),imo(i),iday(i),ihr(i),imin(i),sec,
     &                 rmk1,rmk2,cphase,ipwt,amx,prx,
     &                 xlat,alon,emag(i),depth,e(1,i),
     &                 i1,ifixsolution,i3,eventtype,itrial)
      endif
      if(nobsread.eq.-1)then
         knobs(i)=-1       ! means: end of input-file detected!
         RETURN
      endif
      if(.NOT.(ised.eq.0.or.ised.eq.1.or.ised.eq.2))
     & stop'INPUTDATA: ised flag has no supported value!'
cek
cek check for events with less than three observations
cek
      if(isingle.ne.0) then
        if(nobsread.lt.3) then
            write(6,'(1x,''Event #'',i4)') i
            write(6,*) ' skipped because it has fewer than 3 obs.'
           goto 10001
        endif
      endif
      if(isingle.ne.0)then
         write(2,'(''1 E V E N T   N R .   '',i6,
     &             ''                 '',
     &             ''           0                    0'')')
     &             isingle
         write(6,'(1x,''Event #'',i6)')
     &            isingle
      endif
c
      if(i.le.neqs)then
cEK Dez94  next statement put in effect (ek: i1 set to zero, see above)
         ifx(i)=i1
cEK Dez94 next statement put out of use:
cEK         ifx(i)=0   ! NO ROTATION; do not set dtdr(2) to 0.0  UK87
      else
         icc=i1
      endif
c
c     if ITRIAL = 1 take first station as trial epicenter and
c                   ztrial as trial depth
c     Do the same, if no hypocenter has been read from input !!!
c
      if(itrial.gt.0.or.
     &      (xlat.eq.0.0.and.alon.eq.0.0))then
c
c        find first station:
c
         if(itrial.gt.0)then
            jjmin1=1
            do jj=1,nobsread
               if(ipwt(jj).lt.5)then
                  do jk=1,nsta
                     if(sta(jj).eq.stn(jk)) goto 22222
                  enddo
                  goto 2222   ! jj-th station not on station-list !!!
22222             jjmin1=jj
                  goto 222
               endif
2222           continue
            enddo
 222        jjmin=jjmin1  ! first station in data which is on station-list!!
            do jj=jjmin1,nobsread
               if(ipwt(jj).lt.5)then
                  if(sec(jj).lt.sec(jjmin1))then
                     jjmin=jj
                  endif
               endif
            enddo
cc          write(16,'('' first station is '',a6)') sta(jjmin)
            jjmin1=0
            do jj=1,nsta
               if(sta(jjmin).eq.stn(jj))then
                  jjmin1=jj
                  goto 11111
               endif
            enddo
         endif
11111    continue    ! jj is the first station
c        If no first station found, take 1st in station-list !! :
cuk         if(smn(jj,i).eq.' ') jj=1  ! is wrong!!
         if(jjmin1.eq.0) jj=1
         if(ifixsolution.ne.9)then
            xlat=xla(jj)+0.001
            alon=xlo(jj)+0.001
            if(ifixsolution.ne.1)then
               depth=ztrial
            else
               if(.not.single_turbo) write(16,*)'DEPTH fixed !!!'
               write(6,*) 'DEPTH fixed !'
            endif
         else
            if(.not.single_turbo) write(16,*)'HYPOCENTER fixed !!!'
            write(6,*)'HYPOCENTER fixed !'
            if(icoordsystem.eq.2.and.alon.gt.0.) alon=-alon ! fixed LAT
         endif                                 ! was given in LON E
      endif
c
c     in case the analyst has fixed the depth to less/equal 0.0, he probably
c     wanted to fix it at the surface... set depth to 3km above sea-level;
c     program VELEST will set it properly to the surface!
c
      if(ifixsolution.ne.0.and.itopo.gt.0.and.depth.le.0.0) depth=-3.0
      if(i.gt.neqs)then
         jshot=i-neqs       ! SHOTS
         map2(jshot)=icc
         z=depth
      else
         z=depth+zshift     ! EQS
      endif
c
c     now do transformation 'LAT/LON --> Xkm/Ykm' for trial hypocenter:
c
      if(icoordsystem.eq.2)then
         call GEOKO(e(2,i),e(3,i),xlat,-alon,-1) ! calc. cart. coord.
         e(2,i)=-e(2,i)
      else
         call SDC(e(2,i),e(3,i),xlat,alon,-1) ! calc. cart. coord.
      endif
      e(4,i)=z
      depthsofinput(i)=z
c
      do 15 j=1,nobsread
c
      if(isingle.eq.0)then
         if(ipwt(j).ge.4)then     ! do not accept phase-weights >= 4  !!!
            goto 15
         endif
      else    
         if(ipwt(j).ge.6)then     ! do not accept phase-weights >= 6  !!!
            goto 15               ! but take 4's & 5's for magnitude-calculation
         endif                    ! and 4's also for residual-computation !!!
      endif
      if(nsp.eq.1.and.cphase(j).eq.'s') goto 15
      if(nsp.eq.1.and.cphase(j).eq.'S') goto 15
      if(nsp.eq.1.and.cphase(j).eq.'-') goto 15 ! skip also s-p phases
c
c     if cphase(j).ne. 's' or 'S', a P-phase is assumed, so it's not necessary
c     to check, whether it's a real P or a reflected P ( 'm' or 'M').
c     But if no reflections should be assumed (e.g. ireflector.eq.0),
c     do NOT use any phases marked with 'm' or 'M' !
      if(cphase(j).ne.'p'.and.cphase(j).ne.'P'.and.ipwt(j).lt.5)then
         if(cphase(j).ne.'s'.and.cphase(j).ne.'S')then
            if(cphase(j).ne.'-')then
               if(cphase(j).ne.'m'.and.cphase(j).ne.'M')then	    
                  if(.not.single_turbo)then
                     write(16,*)'WARNING:'
                     write(16,*)'what phase is this ?  ',
     &		          cphase(j),' ???'
		  endif   
                  write(6,*)'what phase is this ?  ',cphase(j),' ???'
                  write(6,*)
                  if(isingle.gt.0)then
                     write(2,'('' DELETED: '',a6,
     &                      '' unknown phase is: '',a1)')
     &                      sta(j),cphase(j)
                     goto 15
               endif  ! if .ne. M
               else ! here, if phase is a reflected one!
                  if(ireflector.eq.0)then
                     i f (.not.single_turbo) t h e n
                        write(16,*)'WARNING:'
                        write(16,*)'subr. INPUTDATA >>> Phase is : ',
     &		             cphase(j)
                        write(16,*)'but ireflector is: ',ireflector
                        write(16,*)'Phase therefore neglected !!'
                     e n d i f  ! single_turbo
                     write(6,*)'subr. INPUTDATA >>> Phase is : ',
     &		          cphase(j)
                     write(6,*)'but ireflector is: ',ireflector
                     write(6,*)'Phase therefore neglected !!'
                     write(6,*)
                     goto 15
                  endif  ! if ireflector
               endif  ! else
	    endif ! if s-p   
         endif  ! if s
      endif  ! if p
c
      do 16 k=1,nsta
      if(sta(j).eq.stn(k)) goto 17
   16 continue
      if(.not.single_turbo)then
         if(isingle.eq.0) write(16,*)'WARNING:     Event # ',i
         if(isingle.gt.0) write(16,*)'WARNING:     Event # ',isingle
         write(16,'('' WARNING:  Station: >>>'',a6,
     &              ''<<< not found in stationlist!'')') sta(j)
         write(16,*)'Phase therefore skipped'
      endif
cc      write(6,*)'Event # ',i
      write(6,'('' WARNING:  Station: >>>'',a6,
     &           ''<<< not found in stationlist!'')') sta(j)
      write(6,*)'Phase therefore skipped'
      write(6,*)
      if(isingle.gt.0) write(2,'('' DELETED: '',a6,
     &                           '' not on station-list'')') sta(j)
      goto 15
17    continue
      ss1=sqrt( (x(k,1)-e(2,i))**2 + (x(k,2)-e(3,i))**2)
      if(ss1.gt.dmax)then
         if(.not.single_turbo)then
            write(16,*)'WARNING:'
            write(16,'('' epicentral distance:'',f6.1,
     &                 '' > dmax ('',f6.1,'') ==> skipping phase !'')')
     &                 ss1,dmax
         endif
         if(isingle.ne.0)then
            write(6,'('' epicentral distance:'',f6.1,
     &                 '' > dmax ('',f6.1,'') ==> skipping phase !'')')
     &                 ss1,dmax
         write(2,'('' DELETED: '',a6,
     &             '' epicentral-distance too large'')') sta(j)
         endif
         goto 15
      endif
c
c     test for only one p-reading of same station per event
c
      if(ipwt(j).eq.5) goto 888  ! dont test for "weight 5" phases
c      
      do ll=1,j-1
         itest=9
         i f (nsp.eq.1.and.(cphase(j).eq.'s'.or.cphase(j).eq.'S'))then
            continue
         e l s e
         if(k.eq.istm(ll,i).and.ipwt(ll).lt.4)then
            if(cphase(j).eq.'p')itest=0
            if(cphase(j).eq.'P')itest=0
            if(cphase(j).eq.'s')itest=1
            if(cphase(j).eq.'S')itest=1
            if(cphase(j).eq.'m')itest=-1
            if(cphase(j).eq.'M')itest=-1
            if(cphase(j).eq.'-')itest=2
            if(sphase(ll,i).eq.itest)then  ! twice the same phase in this event
               if(isingle.eq.0) ie=i
               if(isingle.gt.0) ie=isingle
               if(.not.single_turbo)then
                 write(16,*)'WARNING:'
                 write(16,*)'PHASETEST: POSSIBLE ERROR in phaselist !!'
                 write(16,'(1x,''---> '',3i2.2,1x,2i2.2)')
     &                    iyr(i),imo(i),iday(i),ihr(i),imin(i)
                 write(16,*)'Event=',ie,' Obs-nr. = ',j,' >>> Station ',
     &                      sta(j),' & Phase = ',cphase(j),
     &                      ' already occured!'
               endif
               write(6,*)'PHASETEST: POSSIBLE ERROR in phaselist !!'
               write(6,'(1x,''---> '',3i2.2,1x,2i2.2)')
     &                  iyr(i),imo(i),iday(i),ihr(i),imin(i)
                 write(6,*)'Event=',ie,' Obs-nr. = ',j,' >>> Station ',
     &                      sta(j),' & Phase = ',cphase(j),
     &                      ' already occured!'
               write(6,*)
               if(isingle.eq.0)then
                  if(.not.single_turbo)then
                     write(16,*)'subr. INPUTDATA >>> program will stop'
                     write(16,*)
                  endif
                  iphaseteststopflag=1
               else
                  if(.not.single_turbo)then
                     write(16,*)'INPUTDATA>>> nevertheless, '
     &                           //'program continues'
                     write(16,*)
                  endif
               endif
c              iphaseteststopflag is SAVED at beginning of subr. inputdata
c             ! forget never that error was detected!
c                 stop'PHASE-TEST: FATAL ERROR in phaselist !'
            endif
         endif
         e n d i f
      enddo
c
c     arrive here, if station on stationlist and phase is accepted by all tests:
c
 888  nobs=nobs+1
      amx(nobs)=amx(j)    
      prx(nobs)=prx(j)    ! added 18.9.91 / uk
      if(isingle.gt.0)then
         if(ised.eq.2)then
            prmk(nobs,1)=rmk1(j)
            prmk(nobs,2)=rmk2(j)
         else
            prmk(nobs,1)=' '
            prmk(nobs,2)=' '
         endif
      endif
      do 66 l=1,3
66    d(nobs,l,i)=x(k,l)
      pt(nobs,i)=sec(j)
      kpwt(nobs,i)=ipwt(j)
      istm(nobs,i)=k
      smn(nobs,i)=sta(j)
c***    kpwt(nobs,i) - p or s weight for observation nobs, event i
c***    istm(nobs,i) - station number for observation
c***    smn(nobs,i) - station name for observation
c**    iphase(nobs,i) gives the model number for the observation
      if(nsp.eq.2) goto 21
c
c
        iphase(nobs,i)=model(1)
c
c***     sphase(nobs,i) is 0. for a p observation
c***                       1. for an s observation
cuk*                      -1. for a reflected P observation
chrm                       2. for a S-P phase
c
      sphase(nobs,i)=0.
      if(cphase(j).eq.sc.or.cphase(j).eq.ss) sphase(nobs,i)=1.0
      if(cphase(j).eq.'m'.or.cphase(j).eq.'M') sphase(nobs,i)=-1.0
      if(cphase(j).eq.'-') sphase(nobs,i)=2.0
c
c *** w(nobs,i) is an observation weighting factor :
c
      w(nobs,i)=1.0/(2**(ipwt(j)*2))
c
      if(cphase(j).eq.sc.or.cphase(j).eq.ss.or.cphase(j).eq.'-')
     & w(nobs,i)=swtfac*w(nobs,i)
      if(ipwt(j).gt.4)then
         w(nobs,i)=0.0    ! necessary for single-event-mode !!!
         nobswithw0=nobswithw0+1
      endif
      wsum=wsum+w(nobs,i)
      goto 15
21    continue    !nsp=2
      if(cphase(j).eq.sc.or.cphase(j).eq.ss) goto 22 ! s-phase
      if(cphase(j).eq.'-') goto 22  ! s-p phase
      sphase(nobs,i)=0.
      if(cphase(j).eq.'m'.or.cphase(j).eq.'M') sphase(nobs,i)=-1.0 ! refl. P
      iphase(nobs,i)=model(1)      
c      iphase(nobs,i)=model(2*k-1)
      w(nobs,i)=1.0/(2**(ipwt(j)*2))
      if(ipwt(j).ge.4)then
         w(nobs,i)=0.0    ! necessary for single-event-mode !!!
         nobswithw0=nobswithw0+1
      endif
      wsum=wsum+w(nobs,i)
      goto 15
22    continue
      sphase(nobs,i)=1.
      if(cphase(j).eq.'-') sphase(nobs,i)=2.0 ! s-p phase
c
c     For s-p phases the model parameter is set to the s-wave model
c
c
      iphase(nobs,i)=model(2)
      w(nobs,i)=swtfac*1.0/(2**(ipwt(j)*2))
      if(ipwt(j).ge.4)then
         w(nobs,i)=0.0    ! necessary for single-event-mode !!!
         nobswithw0=nobswithw0+1
      endif
      wsum=wsum+w(nobs,i)
15    continue
c
c     one event is read-in now :
c
 14   continue
      knobs(i)=nobs
cek list P and S obs separateley:
      if(nsp.eq.2) then
         nobsp=0
         nobss=0
         do j=1,nobs
            if(sphase(j,i).eq.0.0) nobsp=nobsp+1
            if(sphase(j,i).eq.1.0) nobss=nobss+1
         enddo
      endif
      if(.not.single_turbo)then
         if(xlat.lt.0.0)then
            cns='S'
            xxlat=-xlat
         else
            cns='N'
            xxlat=xlat
         endif
         if(alon.lt.0.0)then
            cew='E'
            aalon=-alon
         else
            cew='W'
            aalon=alon
         endif
         if(icoordsystem.eq.2)then
            xxx=-e(2,i)
         else
            xxx=e(2,i)
         endif
         if(i.le.neqs)then
cek print P and S number of obs
           if(nsp.eq.2) then
            write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,
     &               1x,f8.4,a1,1x,f6.2,3f7.2,f5.2,i2,i4,3x,2i5)')
     &      i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &      xxlat,cns,aalon,cew,depth,xxx,(e(j,i),j=3,4),
     &      emag(i),ifx(i),knobs(i)-nobswithw0,nobsp,nobss
           else
            write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,
     &                 1x,f8.4,a1,1x,f6.2,3f7.2,f5.2,i2,i4)')
     &      i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &      xxlat,cns,aalon,cew,depth,xxx,(e(j,i),j=3,4),
     &      emag(i),ifx(i),knobs(i)-nobswithw0
           endif

         endif
         if(i.gt.neqs)then  !shots
            write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,
     &                1x,f8.4,a1,1x,f6.2,3f7.2,f5.2,i2,i4)')
     &      i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &      xxlat,cns,aalon,cew,depth,xxx,(e(j,i),j=3,4),
     &      emag(i),map2(i-neqs),knobs(i)-nobswithw0
         endif
      endif
c
c      No more active (format 1607)
c 1607 format(1x,i3,1x,3i2.2,1x,2i2.2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,
c     &       f6.2,3f7.2,f5.2,i2,i4)
c
c     normalize weights:
c
      if(isingle.gt.0)then
         if( (knobs(i)-nobswithw0) .lt. nvar ) RETURN
         if(wsum.le.0.0) RETURN
      endif
      do j=1,nobs
         w(j,i)=w(j,i)*(nobs-nobswithw0)/wsum
      enddo
c
      if(iphaseteststopflag.eq.1.and.i.eq.(neqs+nshot))then
         stop'INPUTDATA >>> PHASE-TEST: FATAL ERROR in phaselist !'
      endif
      RETURN
c
80    continue
      if(.not.single_turbo)then
         write(16,85)
c---- read past end of data
   85    format(' ***** end of data encountered *****')
      write(16,*)'WARNING:  subr. INPUTDATA >>> end of data encountered'
      endif
      write(6,*)'WARNING:'
      stop'subr. INPUTDATA >>> error: end of data!'
c
  998 call OPENERROR('inputdata','EQ-data-input-file FOR008')
  999 call OPENERROR('inputdata','SHOT-data-input-file FOR009')
      return
c
      end ! of subr. inputdata
c
      subroutine STATISTICSOUT
c
c     output the statistics done in this VELEST run
c
      implicit none
      include 'vel_com.f'
c
      integer i,ii,irefr,isour,itotal,nhit,iavgap,lesseq1
      integer mge6,nobslesseq1,nobsmge6,mag,mge5l,mge5i
      integer irefl,idepi,idepl,ntot,i1,i2,i3,i4
      real depthi,depthl,res0,res1,res2,res3,res4
      real tot
      real  rlen,err,mini,maxi,xmag
      integer magnr(49), nobsnr(49)
      integer depthnri(105), depthnrl(105)
      character cstari*51, cstarl*51
c
      do i=1,105
         depthnri(i)=0
         depthnrl(i)=0
      enddo
      do i=1,49
         magnr(i)=0
         nobsnr(i)=0
      enddo
c
      do i=1,neqs
         do ii=1,nltot-1
            if(e(4,i).ge.h(ii).and.e(4,i).lt.h(ii+1))then
               ihypoclayer(ii)=ihypoclayer(ii)+1
            endif
         enddo
         if(e(4,i).ge.h(nltot)) ihypoclayer(nltot)=ihypoclayer(nltot)+1
      enddo
      rlen=0.0
      do i=1,nltot
         rlen=rlen+refraylen(i)       !  rlen  will be total refr.raylength [km]
      enddo
      do i=1,nltot
chrm  The next if else statement seems to be neccessary to avoid
chrm  floating point exceptions in case of rlen = 0
      if (rlen.gt.0.00001) then
         refraylen(i)=100.*refraylen(i)/rlen
      else 
	 refraylen(i) = 0.0
      endif
         if(hitlay(i,1).ge.1)then
            hitlay(i,2)=hitlay(i,2)/hitlay(i,1)  !  average horiz. km in layer i
            hitlay(i,3)=hitlay(i,3)/hitlay(i,1)  !  average verti. km in layer i
         endif
      enddo                                  ! hitlay(i,1) = nofHITS of layer i
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)'RAY-STATISTICS FOR  L A S T  ITERATION'
      write(16,*)'--------------------------------------'
      write(16,*)
      write(16,*)'NHYP : nr of hypocenters in this layer'
      write(16,*)'NREF : nr of headwaves in this layer'
      write(16,*)'%len : % of "refracted km" in this layer '
      write(16,*)'       with respect to all refracted kilometers'
      write(16,*)'NHIT : nr of rays passed thru this layer'
      write(16,*)'xy-km: average horizontal ray length [km] in layer'
      write(16,*)'z-km : average  vertical  ray length [km] in layer'
      write(16,*)'RFLX : number of reflections at bottom of this layer'
      write(16,*)
      write(16,'(a)')' nlay   top ..... bottom     velocity   '//
     &                'NHYP NREF %len  NHIT xy-km  z-km  RFLX'
      write(16,*)
      irefr=0
      isour=0
      irefl=0
      do i=1,nltot-1
         nhit=NINT(hitlay(i,1))
         write(16,'(2x,i2,3x,f6.2,''...'',f6.2,'' km   '',
     &         f5.2,'' km/s'',2(1x,i4),1x,f5.1,1x,i5,2(1x,f5.1),2x,i4)')
     &              i,h(i),h(i)+thk(i),
     &              v(i),ihypoclayer(i),irefrlayer(i),refraylen(i),
     &              nhit,hitlay(i,2),hitlay(i,3),irefllayer(i)
         irefr=irefr+irefrlayer(i)
         isour=isour+ihypoclayer(i)
         irefl=irefl+irefllayer(i)
      enddo
      nhit=NINT(hitlay(nltot,1))
      write(16,'(2x,i2,3x,f6.2,''...'',5x,''  km   '',
     &      f5.2,'' km/s'',2(1x,i4),1x,f5.1,1x,i5,2(1x,f5.1))')
     &      nltot,h(nltot),
     &      v(nltot),ihypoclayer(nltot),irefrlayer(nltot),
     &      refraylen(nltot), nhit,hitlay(nltot,2),hitlay(nltot,3)
      irefr=irefr+irefrlayer(nltot)
      isour=isour+ihypoclayer(nltot)
      irefl=irefl+irefllayer(nltot)
c
      write(16,*)
      write(16,'('' Total nr of events was '',i4)') isour
      write(16,*)
      write(16,'('' Total nr of refracted rays = '',i5)') irefr
      write(16,'('' Total nr of reflected rays = '',i5)') irefl
      write(16,'('' Total nr of   other   rays = '',i5)') noheadwave
      itotal=noheadwave+irefr+irefl
      write(16,'(''                              ------'')')
      write(16,'('' Total nr of    all    rays = '',i5)') itotal
      write(16,*)
      write(16,*)'Average (absolute) error of the raytracers:'
      err=0.0
      if(noheadwave.gt.0) err=(sterr+direrr)/noheadwave
      write(16,'(1x,'' Straight and direct rays : '',f7.2,'' meters'')')
     &           err
      err=0.0
      if(irefr.gt.0) err=refrerr/irefr
      write(16,'(1x,'' Refracted           rays : '',f7.2,'' meters'')')
     &           err
      err=0.0
      if(irefl.gt.0) err=reflerr/irefl
      write(16,'(1x,'' Reflected           rays : '',f7.2,'' meters'')')
     &           err
      write(16,*)
      write(16,*)
      avhraylen=avhraylen/itotal
      write(16,*)'ALL RAYS TOGETHER:'
      write(16,'(1x,''Average horizontal ray length = '',
     &           f6.1,'' km   (Hypocenter --> Station) '')')
     &           avhraylen
      avvraylen=avvraylen/itotal
      write(16,'(1x,''Average  vertical  ray length = '',
     &           f6.1,'' km   (Deepest ray-point --> Station)'')')
     &           avvraylen
      write(16,*)
c
      write(16,*)
      write(16,*)
      write(16,*)'... and some more STATISTICS '
      write(16,*)'---------------------------- '
      write(16,*)' GAP of final epicenters:'
      write(16,*)
      write(16,'(5(''  Event# -> GAP''))')
      write(16,'(5(5x,i3,'' -> '',i3))') (i,igap(i),i=1,neqs)
      write(16,*)
      mini=361
      maxi=-1
      iavgap=0
      do i=1,neqs
         if(igap(i).gt.maxi) maxi=igap(i)
         if(igap(i).lt.mini) mini=igap(i)
         iavgap=iavgap+igap(i)
      enddo
      iavgap=NINT(float(iavgap)/float(neqs))
      write(16,'(1x,''GAPs were between '',i3,'' and '',i3)')
     &           nint(mini),nint(maxi)
      write(16,'(''      (average GAP was '',i3,'')'')') iavgap
      write(16,*)
c
      write(16,*)
      write(16,*)
      write(16,*)' MAGNITUDES of INPUT-DATA:'
      write(16,*)
      write(16,*)'Magnitude (# of events)   ***average number of obs***'
      write(16,*)
      lesseq1=0
      mge6=0
      nobslesseq1=0
      nobsmge6=0
      do i=1,neqs
         mag=NINT(emag(i)*10.)
         if(mag.le.10)then
            nobslesseq1=nobslesseq1+knobs(i)
            lesseq1=lesseq1+1
         endif
         if(mag.gt.10.and.mag.lt.60)then
            magnr(mag-10)=magnr(mag-10)+1
            nobsnr(mag-10)=nobsnr(mag-10)+knobs(i)
         endif
         if(mag.ge.60)then
            mge6=mge6+1
            nobsmge6=nobsmge6+knobs(i)
         endif
      enddo
      cstari='*********1*********2*********3*********4*********5>'
      if(lesseq1.gt.0) nobslesseq1=
     &                 NINT(float(nobslesseq1)/float(lesseq1))
      if(nobslesseq1.gt.50)nobslesseq1=51
      if(nobslesseq1.gt.0)then
         write(16,'(1x,''MAG<= 1.0 '',''('',i3,'') '',a)')
     &              lesseq1,cstari(1:nobslesseq1)
      else
         write(16,'(1x,''MAG<= 1.0 '',''('',i3,'') '',a)')
     &              lesseq1,' '
      endif
      do i=1,49
         xmag=(10.+i)/10.
         if(magnr(i).gt.0) nobsnr(i)=
     &                     NINT(float(nobsnr(i))/float(magnr(i)))
         if(nobsnr(i).gt.50)nobsnr(i)=51
         if(nobsnr(i).gt.0)then
            write(16,'(1x,''MAG = '',f3.1,'' ('',i3,'') '',a)')
     &                xmag,magnr(i),cstari(1:nobsnr(i))
         else
            write(16,'(1x,''MAG = '',f3.1,'' ('',i3,'') '',a)')
     &                xmag,magnr(i),' '
         endif
      enddo
      if(mge6.gt.0) nobsmge6=NINT(float(nobsmge6)/float(mge6))
      if(nobsmge6.gt.50)nobsmge6=51
      if(nobsmge6.gt.0)then
         write(16,'(1x,''MAG>= 6.0 '',''('',i3,'') '',a)')
     &              mge6,cstari(1:nobsmge6)
      else
         write(16,'(1x,''MAG>= 6.0 '',''('',i3,'') '',a)')
     &              mge6,' '
      endif
      write(16,*)
      write(16,*)
      write(16,*)
c
      write(16,*)
      write(16,*)
      write(16,*)' DEPTHs of INPUT-DATA and of LAST ITERATION:'
      write(16,*)' ==========================================='
      write(16,*)
      write(16,*)'      Depth       # of events  '
      write(16,*)
      lesseq1=0
      mge6=0
      mge5l=0
      mge5i=0
      do i=1,neqs
         idepi=NINT(depthsofinput(i))
         if(idepi.lt.100)then
            depthnri(idepi+5)=depthnri(idepi+5)+1
         endif
         if(idepi.ge.100)then
            mge5i=mge5i+1
         endif
         idepl=NINT(e(4,i))
         if(idepl.lt.100)then
            depthnrl(idepl+5)=depthnrl(idepl+5)+1
         endif
         if(idepl.ge.100)then
            mge5l=mge5l+1
         endif
      enddo
      cstari='.........1.........2.........3.........4.........5>'
      cstarl='*********1*********2*********3*********4*********5>'
      do i=1,105
         depthi=i-5
         depthl=i-5
         if(depthnri(i).gt.50) depthnri(i)=51
         if(depthnrl(i).gt.50) depthnrl(i)=51
         if(depthnri(i).gt.0)then
            write(16,'(1x,''DEPTH ( input ) = '',f4.0,'' km :  '',a)')
     &                depthi,cstari(1:depthnri(i))
         else
            write(16,'(1x,''DEPTH ( input ) = '',f4.0,'' km :  '',a)')
     &                depthi,' '
         endif
         if(depthnrl(i).gt.0)then
            write(16,'(1x,''DEPTH (last_IT) = '',f4.0,'' km :  '',a)')
     &                depthl,cstarl(1:depthnrl(i))
         else
            write(16,'(1x,''DEPTH (last_IT) = '',f4.0,'' km :  '',a)')
     &                depthl,' '
         endif
cc         write(16,*)
      enddo
      if(mge5i.gt.50)mge5i=51
      if(mge5l.gt.50)mge5l=51
      if(mge5i.gt.0)then
         write(16,'(1x,''DEPTH ( input ) > 100. '',''km :  '',a)')
     &              cstari(1:mge5i)
      else
         write(16,'(1x,''DEPTH ( input ) > 100. '',''km :  '',a)')
     &              ' '
      endif
      if(mge5l.gt.0)then
         write(16,'(1x,''DEPTH (last_IT) > 100. '',''km :  '',a)')
     &              cstarl(1:mge5l)
      else
         write(16,'(1x,''DEPTH (last_IT) > 100. '',''km :  '',a)')
     &              ' '
      endif
      write(16,*)
      write(16,*)
c
      write(16,*)'Residuals of the stations according to the azimuth:'
      write(16,*)'(RES  = total average residual at station)'
      write(16,*)'(RES1 = average residual of rays from 1st quadrant)'
      write(16,*)'(RES2 = average residual of rays from 2nd quadrant)'
      write(16,*)'(RES3 = average residual of rays from 3rd quadrant)'
      write(16,*)'(RES4 = average residual of rays from 4th quadrant)'
      write(16,*)
      write(16,'(1x,''Stn#  Stn     RES          RES1        '',
     &           '' RES2         RES3         RES4'')')
      do i=1,nsta
         res0=0.0
         res1=0.0
         res2=0.0
         res3=0.0
         res4=0.0
         if(stnazires(i,2).gt.0) res1=stnazires(i,1)/stnazires(i,2)
         if(stnazires(i,4).gt.0) res2=stnazires(i,3)/stnazires(i,4)
         if(stnazires(i,6).gt.0) res3=stnazires(i,5)/stnazires(i,6)
         if(stnazires(i,8).gt.0) res4=stnazires(i,7)/stnazires(i,8)
         tot=stnazires(i,2)+stnazires(i,4)+stnazires(i,6)+stnazires(i,8)
         if(tot.gt.0.0) res0=(
     &                        stnazires(i,2)*res1
     &                       +stnazires(i,4)*res2
     &                       +stnazires(i,6)*res3
     &                       +stnazires(i,8)*res4) / tot
         if(tot.gt.0.0)then
            ntot=NINT(tot)
            i1=NINT(stnazires(i,2))
            i2=NINT(stnazires(i,4))
            i3=NINT(stnazires(i,6))
            i4=NINT(stnazires(i,8))
            write(16,'(1x,i3,3x,a6,1x,5(f7.2,''('',i4,'')''))')
     &      i,stn(i),res0,ntot,res1,i1,
     &      res2,i2,res3,i3,res4,i4
         else
            write(16,'(1x,i3,3x,a6,1x,''-.-'')') i,stn(i)
         endif
      enddo
      write(16,*)
      write(16,*)
c
      RETURN
      end ! of subr. statisticsout
c
      subroutine STATISLOUT
c
c     print location-output in file *.VEL , compatible for program STATISL
c
cek     implemented print output of studentized residuals (studres)
cek  29.3.95 implemented print output of station name for all phases
c             (incl. S)
c
c       principally studentized res = residual/(sigma * sqrt(1.-diagofhat))
c
c       here we print out studres= residual/sqrt(1.-diag of hat matrix)
cek  281093
c
      implicit none
      include 'vel_com.f'
c
      real sec,xlat,xlon,aar,ofd,tfd
      real erh,studres
      real erx,ery,erz
      real xstn,ystn,xhyp,yhyp,azi,tobs,tcorr
      integer nin,idmin,k,js,jd,no,jav,knobs1,i,in,kk,iazi
      integer iamx,iprisecondcard
      character ctime*20, card*90
      character*1 clay, char1, cns,cew
      character*1 q,qs,qd
c
      character*1 class(4)
      data class/'A','B','C','D'/
c
      if( (knobs(1)-nobswithw0) .lt. nvar .and.iabort.eq.0)then
         iabort=1
         if(.not.single_turbo)then
            write(16,*)'knobs(i)-nobswithw0 < nvar !!!'
            write(16,*)'Event cannot be located!!!'
         endif
         write(6,*)'knobs(i)-nobswithw0 < nvar !!!'
         write(6,*)'Event cannot be located!!!'
      endif
      if(iabort.eq.1)then
         write(2,'('' ERROR: insufficient data to locate the quake!'')')
         RETURN
      endif
c
      if(ifixsolution.eq.0)then
         write(2,'(''0 DATE  ORIGIN   TIME   LAT      LON     DEPTH '',
     &             '' MAG  NO  DM GAP  RMS   ALE D-SPR'')')
      endif
      if(ifixsolution.eq.1)then
         write(2,'(''0 DATE  ORIGIN TIME   LAT       LON     *DEPTH*'',
     &             '' MAG  NO  DM GAP  RMS   ALE D-SPR'')')
      endif
      if(ifixsolution.eq.9)then
         write(2,'(''0 DATE  ORIGIN TIME  *LAT*     *LON*    *DEPTH*'',
     &             '' MAG  NO  DM GAP  RMS   ALE D-SPR'')')
      endif
c
      sec=e(1,1)
      nin=imin(1)
      if(sec.lt.0.)then
         sec=sec+60.
         nin=nin-1
      endif
      if(sec.gt.60.)then
         sec=sec-60.
         nin=nin+1
      endif
      if(nin.lt.0)then     !  U.K. 3.Feb.87
         nin=nin+60
         ihr(1)=ihr(1)-1
      endif
c
c     convert hypocenter into degrees:
c
      if(icoordsystem.eq.2)then
         call GEOKO( -e(2,1), e(3,1), xlat, xlon , 1 )  ! calc. LAT/LON
         xlon=-xlon
      else
         call SDC( e(2,1), e(3,1), xlat, xlon , 1 )  ! calc. LAT/LON
      endif
      if(xlat.lt.0.0)then
         cns='S'
         xlat=-xlat
      else
         cns='N'
      endif
      if(xlon.lt.0.0)then
         cew='E'
         xlon=-xlon
      else
         cew='W'
      endif
c
      idmin=999   ! minimum distance (epicenter --> receiver)
      aar=0.0   ! absolute average weighted residual
      do k=1,knobs(1)
         delta=SQRT(  (e(2,1)-d(k,1,1))**2
     &              + (e(3,1)-d(k,2,1))**2  )
         idelta(k)=NINT(delta)
         if(w(k,1).gt.0.0)then
            if(idelta(k).lt.idmin) idmin=idelta(k)
            aar=aar+ABS(res(k,1))*w(k,1)
         endif
      enddo
      aar=aar/(knobs(1)-nobswithw0)
c
c write first line of summary card:
c
      write(2,'(1x,3i2.2,1x,i2,'':'',i2,'':'',f6.3,1x
     &          f7.4,a1,f8.4,a1,1x,f7.3,2x,
     &          f3.1,2x,i2,1x,i3,1x,i3,f5.2,f6.2,
     &          1x,f5.2)')
     &          iyr(1),imo(1),iday(1),ihr(1),nin,sec,
     &          xlat,cns,xlon,cew,e(4,1),
     &          emag(1),knobs(1)-nobswithw0,idmin,igap(1),rms(1),ale(1),
     &          spread
      write(2,'(''0  ERX  ERY  ERZ Q SQD  ADJ  IN NR  '',
     &          ''AVR   AAR  NM AVXM  SDXM IT'')')
c
c     standard deviations :
c     (square-root of diagonalelements of UNIT covariance matrix
c
      erh=sqrt(s(2)**2 + s(3)**2)
      erx=s(2)
      ery=s(3)
      erz=s(4)
      JS=4
      IF((RMS(1).LT.1.00).AND.(ERH.LE.15.0))JS=3
      IF((RMS(1).LT.0.60).AND.(ERH.LE.8.0).AND.(erz.LE.15.0))JS=2
      IF((RMS(1).LT.0.30).AND.(ERH.LE.2.0).AND.(erz.LE.6.0)) JS=1
      JD=4
      no=knobs(1)-nobswithw0
      OFD=e(4,1)      ! focal depth
      TFD=2.*e(4,1)
      IF(OFD .LT.10.) OFD=10.
      IF(TFD .LT. 30.) TFD=30.
      IF((iGAP(1).LE.180).OR.(NO.GE.4).AND.(iDMIN.LE.100)) JD=3
      IF((iGAP(1).LE.135).OR.(NO.GE.5).AND.(iDMIN.LE.TFD )) JD=2
      IF((iGAP(1).LE. 90).OR.(NO.GE.6).AND.(iDMIN.LE.OFD )) JD=1
      JAV=(JS+JD+1)/2
      Q=class(JAV)
      QS=class(JS)
      QD=class(JD)
      knobs1=0
      do i=1,knobs(1)
         if(kpwt(i,1).ge.5) knobs1=knobs1+1
      enddo
      knobs1=knobs(1)-knobs1 ! # of obs with weight less 5 !!!
      IN=0
c
c write second line of summary card
c
      write(2,'(1x,f5.1,f5.1,f5.1,1x,a1,1x,a1,''/'',a1,f6.2,1x,
     &          i2,i3,2f6.2,
     &          i3,2x,f3.1,2x,f3.1,1x,i3)')
     &                          erx,ery,erz,q,qs,qd,steplen,IN,
     &                          knobs1,avres(1),AAR,
     &                          nmag,xmagnitude,sdxmagnitude,nitt
      if(icoordsystem.eq.2)then
         if(nreg.ge.1000)then
            write(2,'(''0 L+T NR:'',i4,1x,a32,''CH-COORD.:'',
     &       f9.3,'' /'',f9.3, '' KM'')') nreg,regionname,-e(2,1),e(3,1)
         else
            write(2,'(''0 F-E NR:'',i4,1x,a32,''CH-COORD.:'',
     &       f9.3,'' /'',f9.3, '' KM'')') nreg,regionname,-e(2,1),e(3,1)
         endif
      else
         write(2,'(''0 F-E NR:'',i4,1x,a32)') nreg,regionname
      endif
      write(2,'(''0 STN  DIST AZM AIN PRMK HRMN  P-SEC  TPOBS  TPCAL '',
     &          '' -TSCOR  P-RES   P-WT IMP STURES'')')
      write(2,'(''        AMX PRX     SRMK XMAG  S-SEC  TSOBS  TSCAL '',
     &          '' -TSCOR  S-RES   S-WT IMP STURES'')')
c
      k=0
c
c  DO LOOP for EACH OBSERVATION
c
      do kk=1,knobs(1)
      k=k+1
c
c     station# are stored in array  ISTM(iobs,ievent)
c     stn(stn#)   is the station-name (a6)
c     smn(iobs,iev)  is the station-name
c     sphase(nobs,iev) : 0 for P, 1 for S, -1 for reflected P (M-phase)
c     and 2 for an s-p phase
c
      clay=' '
      if(sphase(k,1).eq.0) char1='P'
      if(sphase(k,1).eq.1) char1='S'
      if(sphase(k,1).eq.2) char1='-'
      if(sphase(k,1).eq.-1)then
         char1='P'
         clay='M'
      endif
c        Azimuth (stn --> hypoc) = 57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      xstn=d(k,1,1)
      ystn=d(k,2,1)
      xhyp=e(2,1)
      yhyp=e(3,1)
      azi=57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      if(azi.lt.0) azi=azi+360.
c        Azimuth (hypoc --> stn) = Azimuth (stn --> hypoc) + 180 deg
      azi=azi+180.
      azi=MOD(azi,360.)
      iazi=NINT(azi)
c
      tobs=pt(k,1)-e(1,1)  !  P observed travel time
cek
cek correct time for minuite change as previously done for origin time
cek (see above)
cek this correction added by ek 2.May1996
c
      if(tobs.lt.0.)then
         pt(k,1)=pt(k,1)+60.
      endif
      if(tobs.gt.60.)then
         pt(k,1)=pt(k,1)-60.
      endif
c
      tobs=pt(k,1)-e(1,1)  !  P observed travel time
c
      if(tcalc(k).lt.0.)then
         tcalc(k)=tcalc(k)+60.
      endif
      if(tcalc(k).gt.60.)then
         tcalc(k)=tcalc(k)-60.
      endif
c
c station corrections:
      tcorr=ptcor(istm(k,1))  
      if(nsp.eq.2.and.sphase(k,1).eq.1.0) tcorr=stcor(istm(k,1))
      if(nsp.eq.3.and.sphase(k,1).eq.1.0) tcorr=tcorr*vpvs
      if(sphase(k,1).eq.2) tcorr= 0.0
c
      card=' '
      studres=res(k,1)/sqrt(1.-drm(k,k))
      if(studres.gt.999.) studres=999.999
c
cek print P card (or S-card if no P obs for this station available
c
      iprisecondcard=0
      write(card,'(2X,A6,1x,2i4,I4,1X,3A1,I1,a1,2I2,3F7.3,F7.3,1x,
     &             f7.3,1x,f6.2,1x,f6.4,f7.3)')
     &        smn(k,1),idelta(k),iazi,iain(k),prmk(k,1),char1,prmk(k,2),
     &        kpwt(k,1),
     &        clay,ihr(1),nin,pt(k,1),tobs,tcalc(k),
     &        tcorr,res(k,1),w(k,1),drm(k,k),studres
cek      if(w(k,1).eq.0.0) write(card(61:61),'(''*'')')
      if(kpwt(k,1).eq.5)then
         card(21:80)=' ' ! there is no P arrival!
         if(k.lt.knobs(1))then
            if(smn(k+1,1).eq.smn(k,1))then
               write(card(16:19),'(i4)') iain(k+1)  ! angle of next phase (S) !!
               write(2,'(a)') card !next card will be an S
            endif
         else
            goto 99
         endif
      else
         write(2,'(a)') card
      endif
      iamx=NINT(amx(k))
      card=' '
      if(xmagni(k).ne.-13.)then  ! a magnitude was calculated for this observ.!
         write(card(1:15),'(2x,4x,i5,f4.1)') iamx, prx(k)
         write(card(26:28),'(f3.1)') xmagni(k)
      endif
      if(k.lt.knobs(1))then
chrm?   Must anything changed here for s-p phases?  
c
c  NOW  PRINT  S-PHASE 
c   
         if(smn(k+1,1).eq.smn(k,1))then
            if(sphase(k+1,1).eq.1)then
              iprisecondcard=1
              k=k+1
c
              char1='S'
              clay=' '
              tobs=pt(k,1)-e(1,1)
cek
cek correct time for minuite change as previously done for origin time
cek (see above)
cek this correction added by ek 2.May1996
c
              if(tobs.lt.0.)then
                  pt(k,1)=pt(k,1)+60.
              endif
              if(tobs.gt.60.)then
                  pt(k,1)=pt(k,1)-60.
              endif
c
              tobs=pt(k,1)-e(1,1)  !  S observed travel time
c
              if(tcalc(k).lt.0.)then
                  tcalc(k)=tcalc(k)+60.
              endif
              if(tcalc(k).gt.60.)then
                  tcalc(k)=tcalc(k)-60.
              endif
c
              tcorr=ptcor(istm(k,1))
              if(nsp.eq.2.and.sphase(k,1).eq.1.0) then
		 tcorr=stcor(istm(k,1))
              endif
              if(nsp.eq.3.and.sphase(k,1).eq.1.0) then
		 tcorr=tcorr*vpvs
              endif
cek       NOW write a following card (second phase) for same station
cek  29.3.95
              write(card(3:8),'(a6)') smn(k-1,1) !write station name for S
              write(card(23:27),'(3a1,i1,a1)')
     &                          prmk(k,1),char1,prmk(k,2),kpwt(k,1),clay
              studres=res(k,1)/sqrt(1.-drm(k,k))
              if(studres.gt.999.) studres=999.999
              write(card(32:),'(3F7.3,F7.3,1x,f7.3,1x,F6.2,1x,f6.4,
     &                          f7.3)')
     &                          pt(k,1),tobs,tcalc(k),
     &                          tcorr,res(k,1),w(k,1),drm(k,k),studres
cek              if(w(k,1).eq.0.0) write(card(63:63),'(''*'')')
            endif
         endif
      endif
cek next statements may95
      if(iprisecondcard.eq.1) then
         write(2,'(a)') card
      endif
cek      if(kpwt(k,1).lt.5) write(2,'(a)') card
 2022 continue
      if(k.eq.knobs(1)) goto 99
      enddo
c
  99  continue
c      call DATETIME(ctime)  ! get date&time from system
c      write(2,'(''  $$$ '',2x,''VELEST-Version ETH-11FEB92'',
c     &          '' located at: '',a20)') ctime
      write(2,*)
      write(2,*)
c
      return
      end ! of subr. statislout
c
      subroutine NITTOUTPUT(damp)
c
c     output the results of iteration NITT:
c
      implicit none
      real damp
      include 'vel_com.f'
c
      real avdt,avdx,avdy,avdz,aavdt,aavdx,aavdy,aavdz
      real avelo
      integer i,j2,j1,j,k,k2,j11,j22,kj,jjj,ifl,ksta1,k1,kk1
      integer ksta2
      real cc(ist)
      character*1 reflch
c
c     output hypocenters of this iteration:
c
      avdt=0.0
      avdx=0.0
      avdy=0.0
      avdz=0.0
      aavdt=0.0
      aavdx=0.0
      aavdy=0.0
      aavdz=0.0
      write(16,*)
      write(16,'(4h  eq, 7x, 2hot, 5x, 1hx, 6x, 1hy, 6x, 1hz, 6x,3hrms,
     &           4x,5havres,''   dot     dx     dy     dz'' )')
c
      do 29 i=1,legs
c
c     print constrain-info if necessary:
c
      if(isingle.ne.0)then
         if(isconstrain(1).eq.1)then
            write(16,*)' *** nitt<2 --> ',
     &                 'depth-adjustment := 0.0 for event ',i
         endif
         if(isconstrain(2).eq.1)then
            write(16,*)' *** igap>250 --> ',
     &                 'depth-adjustment := 0.0 for event ',i
         endif
         if(isconstrain(3).eq.1)then
            write(16,*)' *** depth-adjustment constrained for event ',i
         endif
      endif
      if(iconstrain(i).eq.1)then
         write(16,'('' ***** depth constrained for event '',i5)') i
      endif
c
c   output new hypocenter-results got in this iteration:
c
      j2=4*i
      j1=j2-3
      if(i.gt.neqs)then  ! for shots
         j1=3*neqs+i
         j2=j1
      endif
      if(icoordsystem.eq.2)then
         write(16,37) i,e(1,i),-e(2,i),(e(j,i),j=3,4),rms(i),avres(i),
     &             b(j1), -b(j1+1), +b(j1+2), b(j1+3)      ! dot  dx  dy  dz
      else
         write(16,37) i,(e(j,i),j=1,4),rms(i),avres(i),
     &                (b(j),j=j1,j2)              ! dot  dx  dy  dz
      endif
 37   format (1x,i4,3x,6f7.2,2x,4f7.3)
      avdt=avdt+b(j1)
      avdx=avdx+b(j1+1)
      avdy=avdy+b(j1+2)
      avdz=avdz+b(j1+3)
      aavdt=aavdt+ABS(b(j1))
      aavdx=aavdx+ABS(b(j1+1))
      aavdy=aavdy+ABS(b(j1+2))
      aavdz=aavdz+ABS(b(j1+3))
c
  29  continue  ! loop over all events (printout)
c
      if(icoordsystem.eq.2) avdx=-avdx
      avdt=avdt/float(legs)
      avdx=avdx/float(legs)
      avdy=avdy/float(legs)
      avdz=avdz/float(legs)
      aavdt=aavdt/float(legs)
      aavdx=aavdx/float(legs)
      aavdy=aavdy/float(legs)
      aavdz=aavdz/float(legs)
      write(16,*)
      write(16,'('' A V E R A G E   of ADJUSTMENTS :'',19x,4f7.3)')
     &           avdt,avdx,avdy,avdz
      write(16,'('' A V E R A G E   of ABSOLUTE ADJUSTMENTS :'',10x,
     &           4f7.3)') aavdt,aavdx,aavdy,aavdz
      write(16,*)
c
c     all new hypocenters printed out now.
c
      if(damp.ne.1.0)then
         write(16,124) damp
124      format(/,' Step length damping of ',f7.5,' was applied.',/)
      else
         write(16,'(/,''NO step length damping applied'',/)')
      endif
c
c     print velocity adjustments and velocity-model here:
c
      if(scale(6).eq.0.0) goto 44  ! no velocity-adjustments...
      j1=4*neqs+nshot+1
      j2=j1+nltot-1
      k=0
      write(16,*)
      write(16,40)
 40   format(' Velocity adjustments:')
      write(16,41)
 41   format(5x,'vp    dvp      hp   reflector')
c
c     do it for all velocity-models
c
      do i=1,nmod
         write(16,71) i
71       format(1x,'Velocity model',i4)
         k2=laysum(i)
         j11=j1+k2-1
         j22=j11+nplay(i)-1
         kj=1
c
c        do it for all layers in this model ( = i )
c
         do jjj=j11,j22
            if(kj.eq.ireflector)reflch=reflchar
            if(kj.ne.ireflector)reflch=' '
            write(16,63) vp(i,kj),b(jjj),hp(i,kj),reflch
 63         format(1x,2f7.3,2x,f7.3,5x,a1)
            kj=kj+1
         enddo
c
c        calculate and print average velocities of the model i :
c
         ifl=1
         write(16,*)
         write(16,*)'Calculation of average velocity starts at layer #',
     &              ifl
         avelo=0
         do kj=ifl+1,nplay(i)
            avelo=avelo + ( hp(i,kj)-hp(i,kj-1) ) * vp(i,kj-1)
            write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &                ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,kj-1),hp(i,kj),vp(i,kj-1),
     &           avelo/(hp(i,kj)-hp(i,ifl)),hp(i,kj)
         enddo
         write(16,*)
c
         ifl=2
         write(16,*)
         write(16,*)'Calculation of average velocity starts at layer #',
     &              ifl
         avelo=0
         do kj=ifl+1,nplay(i)
            avelo=avelo + ( hp(i,kj)-hp(i,kj-1) ) * vp(i,kj-1)
            write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &                ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &           hp(i,kj-1),hp(i,kj),vp(i,kj-1),
     &           avelo/(hp(i,kj)-hp(i,ifl)),hp(i,kj)
         enddo
         write(16,*)
         write(16,*)
      enddo  ! loop over all models
c
c     each velocity-model and its changes printed out now.
c
c     print p- & s- station-corrections here:
c
 44   if(scale(5).eq.0.0) goto 72 ! no station-correction adjustments...
      write(16,*)
      write(16,'('' Adjusted station corrections:'')')
      write(16,48)
 48   format(2x,' stn  ptcor  dpcor  ')
      k1=4*neqs+nshot+nltot+1
      ksta1=ksta
      if(nsp.eq.2) ksta1=ksta/2
      k2=k1+ksta1-1
c
c     print p-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 49
         if(map1(j).gt.ksta1) goto 49
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)             !  STATION-CORRECTION-
cccc         ptcor(j)=ptcor(j)+cc(j)  !  ADJUSTMENT !!! (p-correction)
 49      continue
      enddo
      write(16,50) (stn(j),ptcor(j),cc(j),j=1,nsta)
 50   format(4(2x,a6,2f7.3))
      if(nsp.ne.2) goto 72
      write(16,'('' Adjusted station corrections:'')')
      write(16,73)
73    format(2x,' stn  stcor  dscor  ')
      k1=4*neqs+nshot+nltot+1+ksta1
      ksta2=ksta-ksta1
      k2=k1+ksta2-1
c
c     print s-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 74
         if(map1(j).gt.ksta2) goto 74
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)              ! STATION-CORRECTION-
ccc         stcor(j)=stcor(j)+cc(j)   ! ADJUSTMENT !!! (s-correction)
74       continue
      enddo
      write(16,50) (stn(j),stcor(j),cc(j),j=1,nsta)
c
72    continue
      write(16,*)
c
      RETURN
c
      end ! of subr. nittoutput
c
      subroutine FINALHYPOCOUT
c
c     output all the final hypocenters, print output simultaneous mode
c     and *cnv format
c
      implicit none
      include 'vel_com.f'
c
      real zero,xlat,xlon,xxx,sec
      integer izero,i,nin,j
      real tt(ist)
      character*1 phzz(ist), cns,cew
      parameter (zero=0.0,izero=0)
      integer year19,sec10,lat1000,lon1000,xout,yout,dep10
      integer mag10,rms100,idmin,aa
      real dist,dmin,xsta,ysta
c
      if(.not.single_turbo) write(16,*)
      if(isingle.ne.0)then
         if(.not.single_turbo) write(16,1111)
         write(6,1111)
      else
         if(.not.single_turbo)then
            write(16,11)
         endif
      endif
      do 12 i=1,legs
c
c     convert hypocenter into degrees:
c
      if(icoordsystem.eq.2)then
         call GEOKO( -e(2,i), e(3,i), xlat, xlon , 1 ) ! calc. LAT/LON
         xlon=-xlon
         xxx=-e(2,i)
      else
         call SDC( e(2,i), e(3,i), xlat, xlon , 1 ) ! calc. LAT/LON
         xxx=e(2,i)
      endif
      if(xlat.lt.0.0)then
         cns='S'
         xlat=-xlat
      else
         cns='N'
      endif
      if(xlon.lt.0.0)then
         cew='E'
         xlon=-xlon
      else
         cew='W'
      endif
c
      sec=e(1,i)
      nin=imin(i)
   23 if(sec.lt.0.) goto 13
   15 if(sec.lt.60.) goto 14
      sec=sec-60.
      nin=nin+1
      goto 15
   13 sec=sec+60.
      nin=nin-1
      goto 23
   11 format(1h ,4x,' date    origin   latitude longitude ',   ! file016
     &     ' depth  mag  no  rms      x      y      z')
 1111 format('  date    origin   latitude longitude',          ! screen (unit=6)
     &     '  depth  mag  no  rms      x      y      z')
c
c     output summary information
c
   14 continue
      if(nin.lt.0)then     !  U.K. 3.Feb.87
         nin=nin+60
         ihr(i)=ihr(i)-1
      endif
c
      call GAPCALC(i)
c
      if(.not.single_turbo)then
         write(16,'(1x,i3,1x,3i2.2,1x,2i2.2,f6.2,
     &              1x,f7.4,a1,1x,f8.4,a1,1x
     &              f6.2,f5.2,i4,f6.3,1x,3f7.2)')
     &             i,iyr(i),imo(i),iday(i),ihr(i),nin,sec,
     &             xlat,cns,xlon,cew,e(4,i),
     &             emag(i),knobs(i)-nobswithw0,rms(i),
     &             xxx,(e(j,i),j=3,4)
      endif
c
      if(isingle.ne.0)then
         write(6,'(1x,3i2.2,1x,2i2.2,f6.2,
     &              1x,f7.4,a1,1x,f8.4,a1,1x,
     &              f6.2,f5.2,i4,f6.3,2f7.2,f6.2)')
     &             iyr(i),imo(i),iday(i),ihr(i),nin,sec,
     &             xlat,cns,xlon,cew,e(4,i),
     &             emag(i),knobs(i)-nobswithw0,rms(i),
     &             xxx,(e(j,i),j=3,4)
      endif
c
      if(ismpout.eq.1)then
         year19 = iyr(i) + 1900
	 sec10 = nint(sec * 10.)
	 lat1000 = nint(xlat * 1000.)
	 lon1000 = nint(xlon * 1000.)
	 xout = nint(xxx)
	 yout = nint(e(3,i))
	 dep10 = nint(e(4,i))
	 if (dep10.lt.0.0) then
	    dep10 = 0.0
         endif
	 mag10 = nint(emag(1) * 10.)
	 if (mag10.lt.0.0) then
	    mag10 = 0.0
         endif
	 rms100 = nint(rms(i) * 100.)
c    
c        Calculating dmin
c
         dmin = 9999
         do aa=1,knobs(i)
	    if (kpwt(i,aa).lt.4) then
	       xsta = -x(istm(aa,i),1)
	       ysta = x(istm(aa,i),2)
               dist = sqrt((xxx-xsta)**2 + (e(3,i) - ysta)**2)
	       if (dist.lt.dmin) then
	          dmin = dist
	       endif
	    endif
	 enddo
	 idmin = nint(dmin)
         write(smpline,'(i4,4i2.2,i3.3,i5.5,a1,i6.6,a1,
     &              i3.3,i2.2,''Ml'',i4.4,2i3,''000000000SEDL'',
     &              2i3,i4.4,i2)')
     &              year19,imo(i),iday(i),ihr(i),nin,sec10,
     &              lat1000,cns,lon1000,cew,
     &              dep10,mag10,nreg,xout,yout,
     &              rms100,igap(i),idmin,knobs(i)-nobswithw0
      endif
c
chrm  If running in single mode, the diagonal elements of the covariance matrix
chrm  will also added to the smpline. In the simultaneous mode the line will be
chrm  written to the file here.
c
cek
cek   this option turned off for single event mode since smp-file only in
cek   SED special format
cek
      if(isingle.eq.0.and.ismpout.eq.1) then
         write(11,'(a80)') smpline
      endif
ccc      if(.not.single_turbo) write(16,*)
      if(isingle.ne.0)then
         if(.not.single_turbo)then
          write(16,'(1x,''Event# '',i3,'' GAP = '',i3)')
     &	  isingle,igap(i)
         endif
         if(icoordsystem.eq.2.and.nreg.ge.1000)then
            if(.not.single_turbo)then
               write(16,'(1x,a32,''   L+T Nr.: '',i4)') 
     &      regionname,nreg
            endif
            write(6,'(1x,a32,''   L+T Nr.: '',i4)') 
     &      regionname,nreg
         else
            if(.not.single_turbo)then
               write(16,'(1x,a32,''   F-E Nr.: '',i4)') 
     &         regionname,nreg
            endif
            write(6,'(1x,a32,''   F-E Nr.: '',i4)') 
     &      regionname,nreg
         endif
      endif
      if(icnvout.eq.0) goto 12  ! do NOT write on file07
c
c---- output final hypocenters and travel times to file07
c
cek     write(7,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,2f7.2)')
cek     &          iyr(i),imo(i),iday(i),ihr(i),nin,sec,
cek     &          xlat,cns,xlon,cew,
cek     &          e(4,i),emag(i)
cek
cek next statement changed by E.Kissling, 21.12.90:
cek          added output of gap and rms of event to converted format
cek
      write(7,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,2f7.2,
     &         4x,i3,5x,f5.2)')
     &         iyr(i),imo(i),iday(i),ihr(i),nin,sec,
     &         xlat,cns,xlon,cew,
     &         e(4,i),emag(i),igap(i),rms(i)
c
      imin(i)=nin
      do 18 j=1,knobs(i)
   18 tt(j)=pt(j,i)-e(1,i)
      E(1,I)=SEC
c
      do j=1,knobs(i)
         phzz(j)='P'
         if(sphase(j,i).eq.1.) phzz(j)='S'
         if(sphase(j,i).eq.-1.) phzz(j)='M'
         if(sphase(j,i).eq.2.) phzz(j)='-'
      enddo
      write(7,19) (smn(j,i),PHZZ(j),kpwt(j,i),tt(j),j=1,knobs(i))
      write(7,*)
   19 format(6(a6,a1,i1,f6.2))
c
  12  continue
      if(.not.single_turbo) write(16,*)
c
      RETURN
      end ! of subr. finalhypocout
c
      subroutine OPENERROR(subr,char)
c
      implicit none
      include 'vel_com.f'
c
      character*(*) subr,char
      write(6,*)'WARNING:'
      write(6,*)'SUBROUTINE :',subr,'    ERROR OPENING FILE: ',char
      if(.not.single_turbo)then
         write(16,*)'WARNING:'
         write(16,*)'SUBROUTINE :',subr,'    ERROR OPENING FILE: ',char
      endif
      stop'Openerror; program VELEST stopped !'
      end ! of subr. openerror
c
      subroutine STINP(itime,stn,sfreq,isper,iscon,isdmp,isamp,scor,ier)
c tabelle fuer stationen im input file
c itime= eventtime in minuten
c
c      ifilt(l)= sfreq
c      iseis(l)= isper
c      sconst(l)=float(iscon)/10.
c      sdampf(l)= float(isdmp)/10.
c      voltgain(l)= isamp
c      cormag(l)= float(iscor)/1000.
c
      implicit none
c      
      integer isper,iscon,isdmp,isamp,ier
      real scor
      integer ifirstcall,ilin,ierr,isyr1,ismo1,isdy1,ishr1,ismin1
      integer isyr2,ismo2,isdy2,ishr2,ismin2,iscor
      character*2 sfreq
chrm      character*4 snam,stn
      character*6 stn
cNEW!!:
      character*6 snam6,stn6
      integer itime,istime1,istime2,juliam
c
      character*80 stlin(600)    ! file STLIST.DAT contains currently 600 lines
      integer maxstlin
      save stlin
      save maxstlin, ifirstcall
c
c     for the first time in this subr. --> read seismo-file into STLIN :
c
      if(ifirstcall.ne.10000001)then
         ifirstcall=10000001
         maxstlin=1
1000     read(10,'(a)',end=2000) stlin(maxstlin)
         maxstlin=maxstlin+1
         if(maxstlin.eq.600) stop'STINP>>> array stlin too small!!'
         goto 1000
2000     maxstlin=maxstlin-1
      endif
c
      stn6=stn
      ilin=0
      ierr=0
1     ilin=ilin+1
      if(ilin.gt.maxstlin) goto 40
      read(stlin(ilin),'(a6)') snam6   !!!!NEW
      if(snam6(1:4).ne.stn6(1:4)) goto 1
20    ilin=ilin+1
      if(ilin.gt.maxstlin) goto 40
      read(stlin(ilin),11)isyr1,ismo1,isdy1,ishr1,ismin1,
     &                    isyr2,ismo2,isdy2,ishr2,ismin2,
     &                    sfreq,isper,iscon,isdmp,isamp,scor
11    format(4x,2(i4,4i2,1x),a2,i2,2i3,i6,f6.2)
      if(isyr1.ne.0) then
         istime1=juliam(isyr1,ismo1,isdy1,ishr1,ismin1)
         istime2=juliam(isyr2,ismo2,isdy2,ishr2,ismin2)
         if(itime.lt.istime1.or.itime.gt.istime2) then
            goto 20
         else
30          ilin=ilin+1
            if(ilin.gt.maxstlin) goto 40
            read(stlin(ilin),11) isyr1
            if(isyr1.ne.0) goto 30
         endif
      endif
      iscor= scor*1000
      return
c
c     if station not found in stationfile with seismometer-parameters:
c
40    continue
      ierr= -6
      return
      end ! of subr. stinp
c
      subroutine INPUTCNV(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  i1,i2,i3,eventtype)
c
c     on Output:
c     ---------
c     EQS:   i1=ifx(i)   
c     SHOTS: i1=icc
c
      implicit none
c
      integer maxobsperevent
      parameter (maxobsperevent=180)
      integer nobs
      character*6 sta(maxobsperevent)
      integer*4 itime
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent)
      integer iunit,iwt(maxobsperevent)
      real amx(maxobsperevent), prx(maxobsperevent),
     &     xlat,xlon
      real xmagni, depth, origtime
      real ttime(maxobsperevent)
      character*1 eventtype
      integer i1,i2,i3
      character*1 cns,cew
c
      character*80 cline
      integer iyr,imo,iday,ihr,imin
      integer j,j1,j2
c
  1   eventtype='L'
      j2=-1
      nobs=0
c
   2  read(iunit,'(a)',end=99) cline
      if(cline.eq.' ') goto 2
cek search for end of *.cnv file, marked by 9999
      if(cline.eq.'9999') goto 999

c
cek   next line adjusted by EK 3.12.90
      read(cline,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,
     1      2x,f5.2)',
     &           err=9999,end=999)
     &           iyr,imo,iday,ihr,imin,origtime,
     &           xlat,cns,xlon,cew,depth,xmagni ! i1 =ifx(i)
                !=======================!   trial hypocenter !!
cek  i1-switch is dummy
      i1=0
c
cek     EQS:   i1=ifx(i)     ! no longer in use !!!
c     SHOTS: i1=icc
c
      if(cns.eq.'S') xlat=-xlat
      if(cew.eq.'E') xlon=-xlon
c
      call TIMECLEAR(iyr,imo,iday,ihr,imin,origtime,itime)
c
      j2=0
   10 j1=j2+1
      j2=j1+5
      read(iunit,'(6(a6,a1,i1,f6.2))',end=999)
     &            (sta(j),cphase(j),iwt(j),ttime(j),j=j1,j2)
      do j=j1,j2
         if(sta(j).eq.' ')then
            if(j.eq.1.and.j1.eq.1) goto 1  ! blank-line read !
            j2=j-1   ! event is completely read-in
            goto 99
         endif
      enddo
      goto 10         ! read next input-line
c
   99 continue
      nobs=j2
      do j=1,nobs
         sec(j)=origtime+ttime(j)
      enddo
c
      RETURN
c
  999 continue
      nobs=-1
      RETURN
c
 9999 continue
      write(6,*)'INPUTCNV>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTCNV >>> error!'
c
      end ! of subr. inputcnv
c
      subroutine INPUTARCVEL(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  i1,i2,i3,eventtype)
c
c  implemented by EK 9.11.93 to read velest-archive type input
c  data for single event location mode , replacing older
c  input routine for SED-format. called by ised=1
c
c     on Output:
c     ---------
c     EQS:   i1=ifx(i)   ! no longer in use !!!
c     SHOTS: i1=icc
c
      implicit none
c
      integer maxobsperevent
      parameter (maxobsperevent=180)
      integer nobs
      character*6 sta(maxobsperevent)
      integer*4 itime
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent)
      integer iunit,iwt(maxobsperevent)
      real amx(maxobsperevent), prx(maxobsperevent),
     &     xlat,xlon
      real xmagni, depth, origtime
      real ttime(maxobsperevent)
      character*1 eventtype
      integer i1,i2,i3
      character*1 cns,cew
c
      character*80 cline
      integer iyr,imo,iday,ihr,imin
      integer j,j1,IOstatus
c
  1   eventtype='L'
      nobs=-1
c  ISOSTAT was added by M. Zhang to solve the warning issue (end of the file)
   2  read(iunit,'(a)',end=999,IOSTAT=IOstatus) cline
      if(cline.eq.' ') goto 2
      if(IOstatus > 0) stop
c
      read(cline,'(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,
     1      2x,f5.2)',
     &           err=9999,end=999)
     &           iyr,imo,iday,ihr,imin,origtime,
     &           xlat,cns,xlon,cew,depth,xmagni !  i1  =ifx(i)
                !=======================!   trial hypocenter !!
cek  i1-switch is dummy
      i1=0
c
cek     EQS:   i1=ifx(i)     ! no longer in use !!!
c     SHOTS: i1=icc
c
      if(cns.eq.'S') xlat=-xlat
      if(cew.eq.'E') xlon=-xlon
c
      call TIMECLEAR(iyr,imo,iday,ihr,imin,origtime,itime)
c
      j1=1
      do j=1,maxobsperevent
         read(iunit,'(2x,a6,2x,a1,3x,i1,3x,f6.2)',end=99)
     &            sta(j),cphase(j),iwt(j),ttime(j)
         if(sta(j).eq.' ')then
            j1=j1-1          ! event is completely read-in
            goto 99
         endif
         j1=j
      enddo
c
   99 continue
c      nobs=j1    !!!! Here is a bug, leading to last pick missing.
                  !!!! changed by M. Zhang, July 2021
      nobs=j1+1   
c
      do j=1,nobs
         sec(j)=origtime+ttime(j)
      enddo
c
      RETURN
c
  999 continue
      nobs=-1
      RETURN
c
 9999 continue
      write(6,*)'INPUTARCVEL>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTARCVEL >>> error!'
c
      end ! of subr. inputarcvel
c
c
      subroutine INPUTSED(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  iswt,ifixsolution,ievnr,eventtype)
c
c     read inputfile (phase-list in *.SED format)
c     if NOBS=-1 ==> end of input-file detected!
c
      implicit none
c
      integer maxobsperevent
      parameter (maxobsperevent=180)
      integer nobs
      character*6 sta(maxobsperevent)
      integer*4 itime,itime1(maxobsperevent)
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent),clayP,clayS
      integer iunit,iwt(maxobsperevent)
      real amx(maxobsperevent), prx(maxobsperevent),
     &     xlat,xlon
      real xmagni, depth, origtime
      integer iswt,ifixsolution,ievnr
      character*1 eventtype
c
      character*80 cline
      integer iyr1(maxobsperevent),imo1(maxobsperevent),
     &        iday1(maxobsperevent),ihr1(maxobsperevent),
     &        kmin1(maxobsperevent)
      integer iyr,imo,iday,ihr,imin
      integer j,j1, jjmin,jjmin1
c
      ifixsolution=0
      nobs=0
      j=0
c
   1  read(iunit,'(a)',end=999) cline
      if(INDEX(cline,'INST').gt.0) goto 1
      if(INDEX(cline,'SED').gt.0) goto 1
      if(INDEX(cline,'BOL').gt.0) goto 1
      if(INDEX(cline,'EVENT').gt.0) goto 1
c
      j=j+1
      if(cline(1:4).eq.' ')then
         read(cline,'(17x,i1,i1,f5.2,2f8.3)') iswt,ifixsolution,
     &                                        depth,xlat,xlon
         nobs=j-1
         goto 99  ! one event finished!!!
      endif
      read(cline,'(a6,a1,a1,a1,i1,a1,5i2,f5.2,7x,
     &            f5.2,2x,a1,i1,4x,f3.0,f3.1)',err=9999)
     &            sta(j),rmk1(j),cphase(j),rmk2(j),iwt(j),clayP,
     &            iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),sec(j),
     &            sec(j+1),clayS,iwt(j+1),amx(j),prx(j)
      if(j.eq.1)then
         read(cline(64:64),'(a1)') eventtype
         read(cline(76:79),'(i4)') ievnr
      endif
      if(iwt(j).eq.9)then
         j=j-1
         goto 1   ! p-weight 9 not accepted !(means: only T(s-p) )
      endif
      if(clayP.eq.'m'.or.clayP.eq.'M') cphase(j)='M'
      if(cline(20:24).eq.' ') iwt(j)=5 ! NO P-arrival read !
ccc   VELEST cannot handle reflected S-phases; therefore skip them!
      if(clayS.ne.' ') cline(32:36)=' '    ! <-- set NO S-data if REFLECTED
      if(cline(32:36).ne.' ')then   ! field with observed S-time in sec
         cphase(j+1)='S'
         sta(j+1)=sta(j)
         rmk1(j+1)=' '
         rmk2(j+1)=' '
         iyr1(j+1)=iyr1(j)
         imo1(j+1)=imo1(j)
         iday1(j+1)=iday1(j)
         ihr1(j+1)=ihr1(j)
         kmin1(j+1)=kmin1(j)
         if(iwt(j).eq.5)then ! NO P arrival here; store amx&prx to S-arrival !
            amx(j+1)=amx(j)
            prx(j+1)=prx(j)
         endif
         j=j+1              ! increment J, because two phases have been added!!!
         goto 1
      else
         goto 1   ! no S-phase-data: elements J+1 will be overwritten
      endif
c
   99 continue
      do j=1,nobs
         if(iwt(j).lt.5)then
            jjmin1=j
            goto 222
         endif
      enddo
 222  j=jjmin1
      call TIMECLEAR(iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),sec(j),
     &               itime1(j))
      itime=itime1(j)
      jjmin=jjmin1
      do j=1,nobs
         if(iwt(j).lt.5)then
            call TIMECLEAR(iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),
     &                     sec(j),itime1(j))
            if(itime1(j).lt.itime)then
               itime=itime1(j)
               jjmin=j
            endif
         endif
      enddo
c
c     itime is the earliest arrival-minute;
c     adjust all seconds according to this (earliest) minute
c
  100 continue
      do j1=1,nobs
         if(iwt(j1).lt.5) sec(j1)=(itime1(j1)-itime)*60. + sec(j1)
      enddo
      origtime=sec(jjmin)
      iyr=iyr1(jjmin)
      imo=imo1(jjmin)
      iday=iday1(jjmin)
      ihr=ihr1(jjmin)
      imin=kmin1(jjmin)
c
      RETURN
c
  999 continue   ! come here, if end-of-input-file detected !
      nobs=-1
      RETURN
c
 9999 continue
      write(6,*)'INPUTSED>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTSED >>> error!'
c
      end ! of subr. inputsed
c

      
c
      subroutine INPUTSED_NEW(iunit,nobs,
     &                  sta,iyr,imo,iday,ihr,imin,sec,
     &                  rmk1,rmk2,cphase,iwt,amx,prx,
     &                  xlat,xlon,xmagni,depth,origtime,
     &                  iswt,ifixsolution,ievnr,eventtype,itrial)
c
c     read inputfile (phase-list in *.SED format)
c     if NOBS=-1 ==> end of input-file detected!
c
      implicit none
c
      integer maxobsperevent
      parameter (maxobsperevent=180)
      integer nobs
      character*6 sta(maxobsperevent)
      integer*4 itime,itime_o,itime1(maxobsperevent)
      real sec(maxobsperevent)
      character*1 rmk1(maxobsperevent),rmk2(maxobsperevent),
     &            cphase(maxobsperevent)
      integer iunit,iwt(maxobsperevent)
      real amx(maxobsperevent), prx(maxobsperevent),
     &     xlat,xlon
      real xmagni, depth, origtime
      integer iswt,ifixsolution,ievnr
      character*1 eventtype
c
      character*80 cline
      integer iyr1(maxobsperevent),imo1(maxobsperevent),
     &        iday1(maxobsperevent),ihr1(maxobsperevent),
     &        kmin1(maxobsperevent)
      integer iyr,imo,iday,ihr,imin
      integer j,j1, jjmin,jjmin1
      integer ev_year,ev_month,ev_day,ev_hour,ev_min
      character*8 phase_id,amx_type
      integer     use_flag,amx_flag
      integer     found
      integer     itrial
c
      ifixsolution=0
      nobs=0
      j=0
c
   1  read(iunit,'(a)',end=999) cline
c
c     The following cards will be skipped
c
      if(INDEX(cline,'SED').gt.0) goto 1
      if(INDEX(cline,'BOL').gt.0) goto 1
      if(INDEX(cline,'EVENT').gt.0) goto 1
c
c     The INST card will be partially read
c
      if(INDEX(cline,'INST').gt.0) then
         read(cline,'(29x,i2)') iswt
         goto 1
      endif
c
cccccccccccccccccccccccccccc
c     Hypon input section  c 
cccccccccccccccccccccccccccc
c      read(line,*,err=900)
c     +		latr,lonr,ztr,y0tr,motr,dytr,h0tr,
c     +		m0tr,s0tr,inst
c latr:	trial latitude
c lonr:	  "   longitude
c ztr:	  "   depth
c y0tr:	  "   origin year
c motr:	  "      "   month
c dytr:	  "      "   day
c h0tr:	  "      "   hour
c m0tr:	  "      "   minute
c s0tr:	  "      "   second
c inst: 0=  free solution
c       1=  fix depth
c       2=  fix lat/lon depth
c       3=  fix origin time
c       4=  fix all
ccccccccccccccccccccccccccccccccccccccccccc
c
      if(INDEX(cline,'TRIAL').gt.0) then
         read(cline,*) xlat,xlon,depth,iyr,imo,iday,ihr,imin,
     &   origtime,ifixsolution 
	 xlon = -xlon
	 iyr = iyr - 1900
         if (ifixsolution.eq.4) then 
	    ifixsolution = 9
	 endif
         call TIMECLEAR(iyr,imo,iday,ihr,imin,origtime,itime_o)  	 
c	 
c        Reading the line after the TRIAL card
c
         read(iunit,'(a)',end=999) cline
         read(cline,'(2x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') 
     &        ev_year,ev_month,ev_day,ev_hour,ev_min 
	 read(cline(55:55),'(a1)') eventtype
	 read(cline(71:78),'(i8)') ievnr   
         goto 1
      endif
c
c     If the SKIP card is reached, the input for one event is finished
c
      if(INDEX(cline,'SKIP').gt.0) then
         nobs=j
         goto 99  ! one event finished!!!
      endif
c
c     If here, a phase will be read
c 
      j=j+1
      read(cline,'(a6,4x,a8,a1,a1,f8.3,i2,f5.2,f9.0,1x,a8,i3)',
     &     err=9999)	   
     &     sta(j),phase_id,rmk1(j),rmk2(j),sec(j),use_flag,
     &     prx(j),amx(j),amx_type,amx_flag
      call casefold(phase_id)
      call casefold(amx_type)
      call casefold(rmk1(j))
      call casefold(rmk2(j))
c
c     Treat the supported phase types
c
      found = 0
      if (phase_id.eq.'P       ') then
         found=1
	 cphase(j) = 'P'
      endif
      if (phase_id.eq.'S       ') then
         found=1
	 cphase(j) = 'S'
      endif
      if (phase_id.eq.'PMP     ') then
         found=1
	 cphase(j) = 'M'
      endif
      if (phase_id.eq.'S-P     ') then
         found=1
	 cphase(j) = '-'
      endif
c
c     If the read phase is not supported, it will be skipped
c
      if (found.eq.0) then
         j=j-1
	 goto1
      endif	 
c
c     Conversion of the weights
c
      if (rmk1(j).eq.'I') then
         iwt(j)=0
      endif	 
      if (rmk1(j).eq.'E') then
         iwt(j)=1
      endif	 
      if (rmk1(j).eq.'Q') then
         iwt(j)=2
      endif
      if (rmk1(j).eq.' ') then  ! Should never happen
         iwt(j)=5
      endif
chrm      if (use_flag.eq.0) then
chrm         j=j-1
chrm	 goto 1
chrm      endif
      if (use_flag.eq.0) then
         iwt(j) = 4
      endif
c
c     Assigning year,month,day,hour,min to each phase
c
      iyr1(j)=ev_year
      imo1(j)=ev_month
      iday1(j)=ev_day
      ihr1(j)=ev_hour
      kmin1(j)=ev_min
c
c     Read the next event
c
      goto 1
c
c     If here, all phases of one event are read in
c
c
   99 continue
c
c     Searching for the first phase with weight less than 5
c
      do j=1,nobs
         if(iwt(j).lt.5)then
            jjmin1=j
            goto 222
         endif
      enddo
 222  j=jjmin1
      call TIMECLEAR(iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),sec(j),
     &               itime1(j))
      itime=itime1(j)
      jjmin=jjmin1
      do j=1,nobs
         if(iwt(j).lt.5)then
            call TIMECLEAR(iyr1(j),imo1(j),iday1(j),ihr1(j),kmin1(j),
     &                     sec(j),itime1(j))
            if(itime1(j).lt.itime)then
               itime=itime1(j)
               jjmin=j
            endif
         endif
      enddo
c
c     itime is the earliest arrival-minute;
c     adjust all seconds according to this (earliest) minute
c
  100 continue
      if (itrial.gt.0) then
         do j1=1,nobs
            if(iwt(j1).lt.5) sec(j1)=(itime1(j1)-itime)*60. + sec(j1)
chrm	    itime1(j1) = itime  ! Activate this, if you want to use itime1(j1)
         enddo
         origtime=sec(jjmin)
         iyr=iyr1(jjmin)
         imo=imo1(jjmin)
         iday=iday1(jjmin)
         ihr=ihr1(jjmin)
         imin=kmin1(jjmin)
      else
         do j1=1,nobs
            if(iwt(j1).lt.5)  then
	       sec(j1)=(itime1(j1)-itime_o)*60. + sec(j1)
	       if (sec(j1).lt.0) then
	          sec(j1) = 0
	       endif
	    endif
chrm	    itime1(j1) = itime  ! Activate this, if you want to use itime1(j1)
         enddo  
      endif
c
      RETURN
c
  999 continue   ! come here, if end-of-input-file detected !
      nobs=-1
      RETURN
c
 9999 continue
      write(6,*)'INPUTSED>>> read-error! input-line is:'
      write(6,'(a)') cline
      stop'subr. INPUTSED >>> error!'
c
      end ! of subr. inputsed_new
c
cek    end of vel_io.f
c
cek    begin of vel_kern.f
c
      subroutine SETUNT(nitt,invertratio,nsinv,icount,
     &                       xythet,stathet,othet,vthet,zthet,scale)
c
c     set scale factors
c
      implicit none
      integer nitt,invertratio,nsinv,icount
      real xythet,stathet,othet,vthet,zthet
      real scale(7)
c
c***  set damping so that velocities are adjusted every
c***  "invertratio" iterations
      icount=mod(nitt,invertratio)
c
c     scale(1) : origin-time
c     scale(2) : x
c     scale(3) : y
c     scale(4) : z
c     scale(5) : stacorr
c     scale(6) : veloc.model
c     scale(7) : (not used)
c
      scale(1)=1.0
      scale(7)=1.0   ! not used in this version!!!
      scale(2)=sqrt(othet/xythet)
      scale(3)=scale(2)
      scale(4)=sqrt(othet/zthet)
      scale(5)=0.0
      if(nsinv.ne.0.and.icount.eq.0) scale(5)=sqrt(othet/stathet)
      scale(6)=0.0
      if(icount.eq.0) scale(6)=sqrt(othet/vthet)
      return
      end ! of subr. setunt
c
      subroutine FIXUNT(b,neqs,nshot,nl,ksta,scale,
     &                  vdamp,itotmodels,inltot,nlfirst)
c
c---- restore solution to proper units
c
      implicit none
      integer neqs,nshot,nl,ksta
      real b(4*neqs+nshot+nl+ksta),scale(7)
      integer nlfirst,itotmodels,inltot
      real vdamp(itotmodels,inltot)
      integer i,j,k,l,m
      i=0
      if(neqs.le.0) goto 10
      do 1 j=1,neqs
      do 1 k=1,4
      i=i+1
    1 b(i)=b(i)*scale(k)
10    continue
      if(nshot.le.0) goto 3
      do 2 j=1,nshot
      i=i+1
    2 b(i)=b(i)*scale(7)
3     if(scale(6).eq.0.) return
      l=0
      m=1
      do 5 j=1,nl
      i=i+1
      l=l+1
      if(l.gt.nlfirst) then
         m=m+1
	 l=1
      endif
chrm      write(6,*)m,'  ',l
5     b(i)=b(i)*scale(6)/vdamp(m,l)
      if(ksta.eq.0) return
      do 4 j=1,ksta
      i=i+1
    4 b(i)=b(i)*scale(5)
      return
      end ! of subr. fixunt
c
      subroutine SETUPMATRIXG(neq,i)
c
c     puts one row S(i) of the data kernel into the symmetric matrix G in
c     order to accumulate the normal equations.
c
      implicit none
      integer neq,i
      include 'vel_com.f'
      integer j,nn,nni,mm,is1,kk1,ng,mg,ksta1,ksta2,k,k2
      integer nl,nf,mf,n,k1,is
c
ccc      real WK(8,2),bwk(1)
c
      do 2 j=1,nvar
   2  s(j)=0.
c
c     No trvdrv's for reseted obs:
c
      if(isingle.ne.0.and.w(i,neq).eq.0.0) goto 34
c
c---- calculate hypocenter indices
c***   neq-- event number, i--observation number
      nn=4*neq-2
      nni=nn-1
c***    nni- space for time term; nn thru mm -- 3 places for dtdr's
      if(neq.gt.neqs) nni=3*neqs+neq
      mm=nn+2
      if(neq.gt.neqs) mm=nni
      s(nni)=1.0*scale(1)      
      if (sphase(i,neq).eq.2.0) then ! new: for s-p phases (hrm 30.3.92)
         s(nni) = 0.0      ! set derivative for orig. time to zero
      endif	 
      is1=0
      if(neq.le.neqs) goto 40
      s(nni)=1.0*scale(7)
      kk1=map2(neq-neqs)
      if(nshfix.eq.1.and.kk1.ne.0) s(nni)=0.0
      if(nshcor.eq.0) goto 40
      if(kk1.le.0.or.kk1.gt.ksta) goto 40
      if(scale(5).eq.0) goto 40
c
      if(nsp.eq.3.and.sphase(i,neq).eq.1.) goto 40
      if(nsp.eq.2.and.sphase(i,neq).eq.1.) goto 43
      ng=4*neqs+nshot+nltot+1
      mg=ng-1+ksta
      ksta1=ksta
      if(kk1.gt.ksta1) goto 40
      is1=ng-1+kk1
      goto 44
43    continue
      ng=4*neqs+nshot+nltot+(ksta/2)+1
      mg=4*neqs+nshot+nltot+ksta
      ksta1=(ksta/2)
      ksta2=ksta-(ksta/2)
      if(kk1.gt.ksta2) goto 40
      is1=ng-1+kk1
44    continue
      s(is1)=1.0*scale(5)
40    continue
      k=0
      if(zadj.eq.0.0) dtdr(3)=0.0
cc      if(isingle.ne.0.and.nitt.lt.3) dtdr(3)=0.0
cc      if(isingle.ne.0.and.jgap.gt.250) dtdr(3)=0.0
      if(neq.gt.neqs) goto 3
c
cEK Dez. 1994  next statement again put in effect:
      if(ifx(neq).eq.1) dtdr(2)=0.0   ! no longer in use  U.K. 2.Oct.87
      do 4 j=nn,mm
      k=k+1
   4  s(j)=dtdr(k)*scale(k+1)
   3  continue
c---- calculate velocity indices
      if(scale(6).eq.0.) goto 34
15    k2=iphase(i,neq)
      nl=nplay(k2)
      nf=4*neqs+nshot+laysum(k2)
      mf=nf+nl-1
      k=0
      do 5 n=nf,mf
      k=k+1
      if(veladj.eq.0.0) dtdv(k)=0.0
chrm    5 s(n)=dtdv(k)*scale(6)
    5 s(n)=dtdv(k)*scale(6)/vdamp(k2,k)
c---- calculate indices for station terms
      if(scale(5).eq.0.0) goto 34
      k1=istm(i,neq)
      ksta1=ksta
      if(nsp.eq.2) ksta1=ksta/2
      if(nsp.eq.2.and.sphase(i,neq).eq.1.0) goto 16
      if(nsp.eq.3.and.sphase(i,neq).eq.1.0) goto 34
18    continue
      ng=4*neqs+nshot+nltot+1
      mg=ng-1+ksta
      if(map1(k1).eq.0) goto 34
      if(map1(k1).gt.ksta1) goto 34
      is=ng-1+map1(k1)
      goto 17
16    ng=4*neqs+nshot+nltot+(ksta/2)+1
      ksta2=ksta-(ksta/2)
      mg=ng-1+ksta2
      if(map1(k1).eq.0) goto 34
      if(map1(k1).gt.ksta2) goto 34
      is=ng-1+map1(k1)
17    continue
      s(is)=1.0*scale(5)
34    continue
cc
cc     store one line of G in GG:
cc
c
c      if(isingle.ne.0)then
c        do i1=1,4
c            gg(i,i1)=s(i1)
c        enddo
cc         if(i.eq.knobs(1))then ! matrix GG is fully calculated
cc            call LSVDF(GG,100,knobs(1),4,BWK,1,0,SV,WK,ier)
cc            call ALESUBR(SV,4,ale(1))
ccccccccc  test          ale(1)=-log10( sv(4) )   ! test test
ccccc            write(6,*) sv(1),sv(2),sv(3),sv(4), ale(1)
cc         endif
cc      endif
c
c     now accumulate normal equations (vector s contains the traveltime-derivs.)
c
      call ACCUNORMEQS(s,nvar,res(i,neq),w(i,neq),G,RHT)
c
      return
      end ! of subr. setupmatrixg
c
      subroutine CHECKSOLUTION(istopflag,better)
c
c     checks, whether the current solution is 'better' than the previous one.
c     'better' means here, that data variance has decreased. A warning is
c     output in case that although the data variance has decreased the RMS
c     didn't.
c     If datvar has decreased, the output-variable DECREASING is set to 1  .
c
      implicit none
      include 'vel_com.f'
c
      integer istopflag
      real datvar,xmsqrs2,varat1,varat2
      integer decreasing
      logical better
c
      save datvar,xmsqrs2
c
      decreasing=0
      better=.true.
c
      if(nitt.eq.0)then
         datvar=davar1
         xmsqrs2=xmsqrs1
      endif
      if(.not.single_turbo)then
         write(16,*)
         if(ibackups.eq.0)then
            write(16,'('' Iteration nr '',i2,'' obtained:'')') nitt
         else
            write(16,'(''(Iteration nr '',i2,'')   BACKUP nr '',i1,
     &                 '' obtained:'')') nitt,ibackups
         endif
         write(16,10) davar1,xmsqrs1,sqrt(xmsqrs1)
 10      format(' DATVAR=',f12.6,' mean sqrd residual= ',f12.6,
     &          '  RMS RESIDUAL= ',f12.6)
      endif
c
      if(isingle.ne.0) call GAPCALC(1)
c
      if(isingle.eq.0) call AVRESISTATIST
c
c***  form ratio of old to new variance
      varat1=datvar/davar1
c
      if(isingle.ne.0)then ! stop calc. for single event if changes get small
         if(nitt.gt.2.and.abs(datvar-davar1).lt.1e-6)then
            istopflag=1
            if(.not.single_turbo)then
               write(16,*)'Changes in datvar < 1e-6   : STOPPING...'
            endif
         endif
      endif
c*** form ratio of old to new mean sqrd residual
      varat2=xmsqrs2/xmsqrs1
c*** test to see if variance is increasing
      if(varat1.ge..99)then
         decreasing=1
         goto 1
      endif
c*** test to see if mean sqrd residual is increasing
      if(varat1.lt..99.and.varat2.ge..99)then
c
c       msqrd residual decreasing, but data variance increased:
c
         if(.not.single_turbo)then
        write(16,'('' *** WARNING: the data variance has increased.'')')
         endif
        decreasing=1
        goto 1
      endif
      if(varat2.ge..99)then
         decreasing=1
         goto 1
      endif
c
   1  if(decreasing.eq.1)then
c
c        variance is decreasing:
c
         datvar=davar1
         xmsqrs2=xmsqrs1
      else
c
c        variance increased so backup
c
         better=.false.
      endif
c
      RETURN
c
      end ! of subr. checksolution
      subroutine TRAVELTIME(i,nobs,iresflag)
c
c     computes the forward problem: all the raytracing is done, traveltimes
c     AND traveltime-derivatives are calculated. Moreover, a whole bunch
c     of statistics is done here.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input parameter:                                            c
c                                                                 c
c     i          = ith event                                      c
c     nobs       = nr of obs for the ith event                    c
c     iresflag   =
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer i,nobs,iresflag
      include 'vel_com.f'
c
      integer ii,nittdone,ifirst,k2,nl,k1,mll,nrp,nrtn
      integer nrpdeep,ir,kk1,j,jzz
      real z,ttt,tkh,r1,r2,ster,direr,refrer,refler
      real takeoff_angle
      real dtddrefl,dtdhrefl,pobs,extrat1,extrat2
      integer do_s
      real rp(3,inrpmax),
     &     ss(inrpmax),x2(inrpmax),y2(inrpmax),z2(inrpmax)
      real v1(inltot),vsq1(inltot)
      real dtdr_s(3),dtdv_s(inltot),ttt_s
      save nITTdone, ifirst
      real xxx(1),yyy(1),zzz(1),velvel(1,1,1) !added by M. Zhang in Nov. 2020
      !to solve a compiling issue caused by gfortran upgrading.
c
c     reset statistics-parameter if iteration is new
c     or if NITT is the same but all rays are shot again (e.g. after
c     model/hypocenter-backup)
c
      if(isingle.gt.0)then
         if(kpwt(nobs,i).eq.5)then  ! this is not a real phase!!!
            if(w(nobs,i).ne.0.0) stop'TRAVELTIME>>> error in obs-wt!'
            res(nobs,i)=0.0
            do ii=1,3
               dtdr(ii)=0.0
            enddo
            RETURN
         endif
      endif
      if(isingle.eq.0)then
         if(i.eq.1.and.nobs.eq.1.and.nitt.gt.nITTdone)ifirst=1 ! 1st obs 1x done
         if(i.eq.1.and.nobs.eq.1.and.nitt.eq.nITTdone)ifirst=2
c ifirst=1   --->  first of all the obs for the first time done in this iterat.
c ifirst=2   --->  first of all the obs for the 2nd time coming in this iterat.
         if(nitt.gt.nITTdone.or.ifirst.eq.2)then
            nITTdone=nitt
            ifirst=1
            call RESETSTATIS
         endif
      endif
c
c     Reset s-p switch
c
      do_s = 1
c
c
c    z -- event depth
      z=e(4,i)
c    k2 -- velocity model number
      k2=iphase(nobs,i) 
   10 continue      
c     nl -- number of layers for model k2
      nl=nplay(k2)
c    k1 -- station number
      k1=istm(nobs,i)
c    v,vsq are the velocity and velocity squared
c     h(ii) is the depth to the top of layer ii
c     thk(ii) -- thickness of layer ii
c
      do ii=1,nl
         v(ii)=vp(k2,ii)
         v1(ii)=vp(k2,ii)
         if(nsp.eq.3.and.sphase(nobs,i).eq.1.0) v(ii)=vp(k2,ii)/vpvs
         if(nsp.eq.3.and.sphase(nobs,i).eq.2.0) then
	    if (do_s.eq.1) then
	       v(ii)=vp(k2,ii)/vpvs 
	    endif
	 endif   
         vsq(ii)=v(ii)**2
         vsq1(ii)=v1(ii)**2
         h(ii)=hp(k2,ii)
         thk(ii)=thkp(k2,ii)
      enddo
c
c     calculate tkh, r1,r2,delta :
c
      ttt=0.
      tkh=x(k1,3)-h(1)   ! distance from station up to top of model
      r1=(e(2,i)-x(k1,1))
      r2=(e(3,i)-x(k1,2))
c check if event too shallow
      if(itopo.gt.0)then
         if(e(4,i).lt.0.0)then  ! depth above sea-level...so near surface !
            call CHTOP(-e(2,i),e(3,i),zmin,
     &                 topo1file,topo2file) ! zmin:==surface at this point
         else
            zmin=zmininput ! depth below zero ... so zmin=zmininput
         endif
      endif
      if(i.le.neqs)then
         if(e(4,i).le.zmin)then
            if(isingle.eq.0) e(4,i)=zmin+0.1     ! earthquakes
            if(isingle.ne.0) e(4,i)=zmin+0.011   ! earthquake, single event
         endif
      else
         if(e(4,i).le.zmin) e(4,i)=zmin+0.011   ! shots
      endif
      if(ifixsolution.gt.0.and.e(4,i).le.0.0) e(4,i)=zmin+0.001 ! fix depth to
c                                                        !   min_depth allowed!!
cc      r3=(e(4,i)-x(k1,3))
      delta=sqrt(r1*r1+r2*r2)
c
c     set reflector-layer-boundary :
c
      MLL=0
      if(ireflector.gt.0)then
         if(sphase(nobs,i).eq.-1.0)then
            MLL=ireflector    ! P is reflected at bottom of layer ireflector
         endif
      endif
c
c     do the ray-tracing :
c
      xxx(1)=1.0
      yyy(1)=1.0
      zzz(1)=1.0
      velvel(1,1,1)=1.0
c      call RAYPATH(1,1,1,1.,1.,1.,1.,nl,thk,h,v,vsq, ! old one
      call RAYPATH(1,1,1,xxx,yyy,zzz,velvel,nl,thk,h,v,vsq,
     &            e(2,i),e(3,i),e(4,i),x(k1,1),x(k1,2),x(k1,3),
     &            rp,nrp,nrtn,jl,tkj,1,ttt,MLL,ster,direr,refrer,refler,
     &            DTDDrefl,DTDHrefl)
c
      call CHECKRAYPATH(rp,nrp)
c
c
c     raypoints stored in RP [1st index: x,y,z; 2nd: raypt# (max. 2*nlay)]
c     nrp = nr of raypoints stored at moment in RP
c
      if(itopo.eq.2)then !check each raypoint whether it is below surface or not
        call RAYPOINTCHECK(rp,nrp,stn(k1))
      endif
c
c     if ray travels through the air, bend it below surface !!!
c
      if(itopo.eq.3) call BENDRAY(rp,nrp,stn(k1),v(1),ttt)
c
c     sum raytracer-errors
c
c     
c     For s-p phases the next statements should be done after the
c     p-wave calculation
c
      if (sphase(nobs,i).ne.2.0.or.do_s.eq.0) then     
         if(isingle.eq.0)then
            sterr=sterr+ster
            direrr=direrr+direr
            refrerr=refrerr+refrer
            reflerr=reflerr+refler
         endif
c      
         if(isingle.eq.0) call LAYERHIT(rp,nrpdeep,nl,nrp,mll)
c      
         if(irayout.eq.1)then
            write(13,1301) i,stn(k1),nrp
            if(icoordsystem.eq.2)then
               write(13,1302) (-rp(1,ir),rp(2,ir),rp(3,ir),ir=1,nrp)
            else
               write(13,1302) (rp(1,ir),rp(2,ir),rp(3,ir),ir=1,nrp)
            endif
 1301    format(/,' event =',i5,' station ',a6,' num. raypoints= ',i5)
 1302    format(3(2x,3f7.2))
         endif
      endif ! if .NOT. ...	 
c
c     calculate traveltime-derivatives according to the raytype :
c
      goto (5,5,6,6,5,6,5,8),nrtn   !   nrtn=8  if ray is reflected !!
c
      stop'TRAVELTIME>>> illegal nrtn from raytracer!!!'
c
  5   call TRAVDERIV('direct',
     &                nl,mll,v1,vsq1,rp,nrp,x2,y2,z2,ss,r1,r2,i,nobs)
      goto 40
  6   call TRAVDERIV('refracted',
     &                nl,mll,v1,vsq1,rp,nrp,x2,y2,z2,ss,r1,r2,i,nobs)
      goto 40
  8   call TRAVDERIV('reflected',
     &                nl,mll,v1,vsq1,rp,nrp,x2,y2,z2,ss,r1,r2,i,nobs)
 40   continue
c
c     For s-p phases, the traveltime of the s-wave was calculated. Now
c     the p-wave calculation will be performed.
c
      if (sphase(nobs,i).eq.2.0.and.do_s.eq.1) then
         do_s = 0
	 ttt_s = ttt
	 do ii=1,3
	   dtdr_s(ii) = dtdr(ii)
	 enddo  
	 do ii=1,nl
	    dtdv_s(ii) = dtdv(ii)
	 enddo   
	 if(nsp.eq.2) then
	    k2 = k2 - 1
	 endif   
	 goto 10
      endif
      if(sphase(nobs,i).eq.2.0) then
c
c        calculate the s-p traveltime derivative 
c
         do ii=1,3
            dtdr(ii)=dtdr_s(ii) - dtdr(ii)
         enddo
c
c        calculate the s-p time difference
c
         ttt = ttt_s - ttt
      endif	 
c
c     calculate residual by using appropriate station-corrections :
c
      if (sphase(nobs,i).ne.2.0) then
         pobs=pt(nobs,i)-e(1,i)
      else
         pobs=pt(nobs,i)
      endif	 
      extrat1=ptcor(k1)
      extrat2=0.0
      if(nsp.eq.2.and.sphase(nobs,i).eq.1.0)
     &extrat1=stcor(k1)
      if(nsp.eq.2.and.sphase(nobs,i).eq.2.0) then  ! s-p phases
         extrat1=stcor(k1) - ptcor(k1)
      endif	 
      if(nsp.eq.3.and.sphase(nobs,i).eq.1.0)
     &extrat1=ptcor(k1)*vpvs
      if(nsp.eq.3.and.sphase(nobs,i).eq.2.0) then ! s-p phases
         extrat1=ptcor(k1)*vpvs - ptcor(k1)
      endif
      if(nshcor.eq.0) goto 600
      if(i.le.neqs) goto 600
      kk1=map2(i-neqs)
      if(kk1.eq.0) goto 600
      do 601 j=1,nsta
      if(kk1.eq.map1(j)) goto 602
601   continue
      goto 600
c
602   extrat2=ptcor(j)
      if(nsp.eq.2.and.sphase(nobs,i).eq.1.0)
     &extrat2=stcor(j)
      if(nsp.eq.3.and.sphase(nobs,i).eq.1.0)
     &extrat2=ptcor(j)*vpvs
600   continue
      res(nobs,i)=pobs-(ttt+extrat1+extrat2)
c
      if(isingle.ne.0)then
         tcalc(nobs)=ttt
         takeoff_angle=
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2 )
     &   /
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2
     &        +(rp(3,2)-rp(3,1))**2 )
         takeoff_angle=57.296*ASIN(takeoff_angle)
         if( (rp(3,2)-rp(3,1)) .lt. 0.0 )then
            takeoff_angle=180.-takeoff_angle     ! ray is going upwards
         endif
         iain(nobs)=NINT(takeoff_angle)  ! with respect to positive z; downwards
      endif
c
c---- save residual according to the ray-type:
c
      if(isingle.eq.0) call RESISAVE(nrp,nrpdeep,rp,nobs,i,k1,mll)
c
c
c     if VELEST is used in single-event-mode (isingle <> 0 ),
c     set weight to 0 if abs(residual) after 2nd iteration is still > 2. sec
c     ... and normalize the weights once more !
c     --> if a reading with a reseted weight gets a residual which has
c         become small again, 'revive' it !
c
cek      if(isingle.ne.0)then
cek         if(nitt.gt.2.and.abs(res(nobs,i)).gt.2.0
cek     &               .and.igap(i).lt.250)then
cek            call REJECTOBS(i,nobs,iresflag)
cek            call GAPCALC(1)
cek         endif
cek         if(nitt.gt.2.and.abs(res(nobs,i)).lt.1.0
cek     &               .and.igap(i).lt.250)then
cek            if(w(nobs,i).eq.0.0)then
cek               call REVIVEOBS(i,nobs,iresflag)
cek               call GAPCALC(1)
cek            endif
cek         endif
cek      endif
c
      tctime(nobs,i)=ttt
      h(1)=x(k1,3)-tkh ! reset h(1) after it has been altered in subr. raypath!!
      thk(1)=thk(1)+tkh
      if(idrvout.eq.1)then
         write(21,*)'i=',i,'  nobs=',nobs,'  nitt=',nitt,'  nrtn=',nrtn
         write(21,500) ttt,pobs,res(nobs,i)
         write(21,501) (dtdr(jzz),jzz=1,3)
         write(21,502) (dtdv(jzz),jzz=1,nl)
 500     format(' ttt=',f10.5,' pobs= ',f10.5,' res= ',f10.5)
 501     format('dtdr=',3f10.5)
 502     format('dtdv=',7f10.5)
      endif
510   continue
      return
      end ! of subr. traveltime
c
      subroutine RESETSTATIS
c
c     called by subr. TRAVELTIME. The statistics-variables are reset here.
c
      implicit none
      include 'vel_com.f'
      integer ihitl,jhitl
c
      if(irflout.eq.1) rewind(77)
      if(irfrout.eq.1) rewind(78)
      if(iresout.eq.1) rewind(79)
c
      do ihitl=1,nltot
         irefllayer(ihitl)=0
         irefrlayer(ihitl)=0
         refraylen(ihitl)=0.0
         hitlay(ihitl,1)=0.0
         hitlay(ihitl,2)=0.0
         hitlay(ihitl,3)=0.0
      enddo
      noheadwave=0     ! nr of straight & direct waves
      avhraylen=0.0    ! average horiz. raylength
      avvraylen=0.0    !    "    vertic.    "
      sterr=0.0        ! raytracer-errors
      direrr=0.0       !    "
      refrerr=0.0      !    "
      reflerr=0.0      !    "
      avrefrres=0.0    ! average residual for this ray-type
      avotheres=0.0    !    "
      avreflres=0.0    !    "
      abrefrres=0.0    ! average absolute residual for this ray type
      abotheres=0.0    !    "
      abreflres=0.0    !    "
      nrrefrres=0.0    ! nr of residuals for this ray-type
      nrotheres=0.0    !    "
      nrreflres=0.0    !    "
      do ihitl=1,nsta
         do jhitl=1,8
            stnazires(ihitl,jhitl)=0.0
         enddo
      enddo
c
      RETURN
      end ! of subr. resetstatis
c
      subroutine RESISAVE(nrp,nrpdeep,rp,nobs,i,k1,mll)
c
c     saves the residual according to the ray-type
c
      implicit none
      integer nrp,nrpdeep,nobs,i,k1,mll
      include 'vel_com.f'
c
      real xxx,yyy,xhyp,yhyp,xstn,ystn,azi,dist
      integer iazi
      real rp(3,inrpmax)
c
      if(rp(3,nrpdeep).eq.rp(3,nrpdeep+1))then ! ray is horizontal : refracted !
         avrefrres=avrefrres+    res(nobs,i)       ! refracted
         abrefrres=abrefrres+ABS(res(nobs,i))      ! refracted
         nrrefrres=nrrefrres+1
c
         if(lmax.eq.10)then ! refractor is Moho...
            if(irfrout.eq.1)then
c
c              write (x,y) of point A, where ray enters MOHO,
c              of point B, where ray leaves MOHO and of point between A and B
c
               xxx=rp(1,nrpdeep)
               if(icoordsystem.eq.2) xxx=-xxx
               yyy=rp(2,nrpdeep)
               write(78,'(2x,''RFR'',3f10.3,''        .1'')')
     &                    xxx,yyy,res(nobs,i)
               xxx=(rp(1,nrpdeep)+rp(1,nrpdeep+1))/2.
               if(icoordsystem.eq.2) xxx=-xxx
               yyy=(rp(2,nrpdeep)+rp(2,nrpdeep+1))/2.
               write(78,'(2x,''RFR'',3f10.3,''        .1'')')
     &                    xxx,yyy,res(nobs,i)
               xxx=rp(1,nrpdeep+1)
               if(icoordsystem.eq.2) xxx=-xxx
               yyy=rp(2,nrpdeep+1)
               write(78,'(2x,''RFR'',3f10.3,''        .1'')')
     &                    xxx,yyy,res(nobs,i)
             endif
         endif
c
      else
         if(mll.eq.0)then
            avotheres=avotheres+    res(nobs,i)    ! straight or direct
            abotheres=abotheres+ABS(res(nobs,i))   ! straight or direct
            nrotheres=nrotheres+1
         else
            avreflres=avreflres+    res(nobs,i)    ! reflected
            abreflres=abreflres+ABS(res(nobs,i))   ! reflected
            nrreflres=nrreflres+1
c
            if(irflout.eq.1)then
              if(icoordsystem.eq.2)then
                 xxx=-rp(1,nrpdeep)
              else
                 xxx=rp(1,nrpdeep)
              endif
              write(77,'(2x,''RFL'',3f10.4,''        .1'')')
     &                   -rp(1,nrpdeep),rp(2,nrpdeep),res(nobs,i)
            endif
         endif
      endif
c
c     save residual for station# k1 according to the azimuth:
c
      xhyp=rp(1,1)
      yhyp=rp(2,1)
      xstn=rp(1,nrp)
      ystn=rp(2,nrp)
c        Azimuth (stn --> hypoc) = 57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      azi=57.296*ATAN2(xhyp-xstn,yhyp-ystn)
      if(azi.lt.0) azi=azi+360.
      iazi=0
      if(azi.ge.  0.0 .and. azi.lt. 90.0) iazi=1
      if(azi.ge. 90.0 .and. azi.lt.180.0) iazi=2
      if(azi.ge.180.0 .and. azi.lt.270.0) iazi=3
      if(azi.ge.270.0 .and. azi.le.360.0) iazi=4
      stnazires(k1,2*iazi-1)=stnazires(k1,2*iazi-1)+res(nobs,i)
      stnazires(k1,2*iazi)=stnazires(k1,2*iazi)+1.
c
c     type residuals as a function of focus-receiver distance:
      if(iresout.eq.1)then
         dist=(  (rp(1,1)-rp(1,nrp))**2
     &          +(rp(2,1)-rp(2,nrp))**2
     &          +(rp(3,1)-rp(3,nrp))**2  )
         dist=sqrt(dist)
         write(79,'(1x,f6.2,2x,f7.3)') dist,res(nobs,i)
      endif
c
      RETURN
      end ! of subr. resisave
c
      subroutine ADJUSTMODEL(damp)
c
c     adjust model-vector by the solution just obtained
c
      implicit none
      real damp
      include 'vel_com.f'
c
      integer jj,i,n,iminold,j,itime,j1,j2,k,k2,j11,j22
      integer juliam
      integer kj,jjj,k1,ksta1,kk1,ksta2
      real cc(ist)
      real dmin
c
c     adjust hypocenters of this iteration and output them:
c
      jj=0
c
      do 29 i=1,legs
      n=4
      if(i.gt.neqs) n=1   ! for shots adjust only the origin-time !!!
ccc      call TIMECLEAR(iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),iminold)
      if(iyr(i).lt.100)then
         iminold=JULIAM(iyr(i)+1900,imo(i),iday(i),ihr(i),imin(i))
      else
         iminold=JULIAM(iyr(i),imo(i),iday(i),ihr(i),imin(i))
      endif
      do j=1,n
         jj=jj+1
         e(j,i)=e(j,i) + b(jj)
      enddo
      if(e(1,i).lt.0.0)then   ! change reference minute
        call TIMECLEAR(iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i),
     &                                                      itime)
        do j=1,knobs(i)
           pt(j,i)=pt(j,i)+(iminold-itime)*60.
        enddo
      endif
      do j=1,3
         isconstrain(j)=0
      enddo
      iconstrain(i)=0
c
c     constrain focal depth if necessary
c
      if(isingle.ne.0)then
         if(nitt.lt.2)then
            isconstrain(1)=1
c            write(16,*)' *** nitt<2 --> ',
c     &                 'depth-adjustment := 0.0 for event ',i
            e(4,i)=e(4,i)-b(4)
            b(4)=0.0
         endif     
         if(igap(1).gt.250)then
            dmin=999.9   ! minimum distance (epicenter --> receiver)
            do k=1,knobs(1)
               delta=SQRT(  (e(2,1)-d(k,1,1))**2
     &                    + (e(3,1)-d(k,2,1))**2  )
               if(delta.lt.dmin) then
	          dmin=delta
	       endif
            enddo
	    if (dmin.gt.15.) then
               isconstrain(2)=1
               if(iconstrain(i).eq.1) iconstrain(i)=3
c              write(16,*)' *** igap>250 --> ',
c     &                  'depth-adjustment := 0.0 for event ',i
               e(4,i)=e(4,i)-b(4)
               b(4)=0.0
	    endif
         endif
         if(abs(b(4)).gt.zadj)then
            isconstrain(3)=1
c            write(16,*)' *** depth-adjustment constrained for event ',i
            if(b(4).gt.0.0)then
               e(4,i)=e(4,i)-b(4)+zadj
               b(4)=zadj
            else
               e(4,i)=e(4,i)-b(4)-zadj
               b(4)=-zadj
            endif
         endif
      else    ! simultanous inversion
         if(i.le.neqs.and.abs(b(jj)).gt.zadj)then ! b(jj)=depth-adj. if event i
            iconstrain(i)=1
c            write(16,*)' *** depth-adjustment constrained for event ',i
            if(b(jj).gt.0.0)then
               e(4,i)=e(4,i)-b(jj)+zadj
               b(jj)=zadj
            else
               e(4,i)=e(4,i)-b(jj)-zadj
               b(jj)=-zadj
            endif
         endif
      endif
c      effdeltaz(i)=0.0
      if(itopo.gt.0)then
         if(e(4,i).lt.0.0)then  ! depth above sea-level...so near surface !
            call CHTOP(-e(2,i),e(3,i),zmin,
     &                 topo1file,topo2file) ! zmin:==surface at this point
         else
            zmin=0.0 ! depth below zero ... so zmin:==0.0 is sufficiant !
         endif
      endif
      if(e(4,i).lt.zmin.or.(ifixsolution.gt.0.and.e(4,i).le.0.0))then
c         effdeltaz(i)=b(jj)-(e(4,i)-zmin)
c
c        instead of calculating effdeltaz constrain adjustment-vector-element
c        b(jj) directly:
c
         b(jj)=b(jj)-(e(4,i)-zmin)
c
         e(4,i)=zmin
         iconstrain(i)=1
c         write(16,'('' ***** depth constrained for event '',i5)') i
      endif
c
  29  continue  ! loop over all events (adjust)
c
c     all new hypocenters adjusted now.
c
c
c  do velocity adjustments here:
c
      if(scale(6).eq.0.0) goto 44  ! no velocity-adjustments...
      j1=4*neqs+nshot+1
      j2=j1+nltot-1
      k=0
      if(.not.single_turbo)then
         write(16,*)
         write(16,40)
 40      format(' doing velocity adjustments now...')
      endif
c
c     do it for all velocity-models
c
      do 4242 i=1,nmod
      k2=laysum(i)
      j11=j1+k2-1
      j22=j11+nplay(i)-1
      kj=1
c
c     do it for all layers in this model ( = i )
c
      do 42 jjj=j11,j22
      vp(i,kj)=vp(i,kj)+b(jjj)   !   <<---- VELOCITY-ADJUSTMENT
c
c     If reflected phases are used, the velocity of the layer just above
c     the reflector is NOT allowed to have a lower velocity than any of
c     the above layers!
c     Because usually the reflector is almost at the bottom of the model,
c     namely the MOHO, we simplify the condition and don't allow any
c     low-velocity-layers at all in the whole model (of course  o n l y
c     if reflected phases are used in this VELEST-run !)
      if(lowveloclay.eq.0)then
         if(kj.gt.1)then
            if(vp(i,kj).lt.vp(i,kj-1))then
               i f (.not.single_turbo)then
               write(16,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') kj
               write(16,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               endif
               write(6,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') kj
               write(6,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               write(6,*)
               b(jjj)=0.0
               vp(i,kj)=vp(i,kj-1)+0.001
            endif
         endif
      endif
c
      kj=kj+1
42    continue
c
4242  continue
c
c     each velocity-model adjusted now.
c
c     adjust p- & s- station-corrections here:
c
 44   if(scale(5).eq.0.0) goto 72 ! no station-correction adjustments...
      if(.not.single_turbo)then
         write(16,*)
         write(16,'('' doing station-correction adjustments...'')')
         write(16,*)
      endif
      k1=4*neqs+nshot+nltot+1
      ksta1=ksta
      if(nsp.eq.2) ksta1=ksta/2
      k2=k1+ksta1-1
c
c     do p-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 49
         if(map1(j).gt.ksta1) goto 49
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)             !  STATION-CORRECTION-
         ptcor(j)=ptcor(j)+cc(j)  !  ADJUSTMENT !!! (p-correction)
 49      continue
      enddo
      if(nsp.ne.2) goto 72
      k1=4*neqs+nshot+nltot+1+ksta1
      ksta2=ksta-ksta1
      k2=k1+ksta2-1
c
c     do s-correction for all stations:
c
      do j=1,nsta
         cc(j)=0.
         if(map1(j).eq.0) goto 74
         if(map1(j).gt.ksta2) goto 74
         kk1=k1-1+map1(j)
         cc(j)=b(kk1)              ! STATION-CORRECTION-
         stcor(j)=stcor(j)+cc(j)   ! ADJUSTMENT !!! (s-correction)
74       continue
      enddo
c
72    continue
c
      RETURN
c
      end ! of subr. adjustmodel
c
      subroutine RMSDATVAR
c
c     compute RMS and DATVAR for all the events (EQs & shots)
c
      implicit none
      include 'vel_com.f'
      integer jj,i,j2,j,knobst
      real tres
c
c
c   now calculate rms for all events
c
      jj=0
      davar1=0.0
      xmsqrs1=0.0
      tres=0.0
      knobst=0
      do 11 i=1,legs
         j2=knobs(i)
         knobst=knobst+j2
         do 12 j=1,j2
            if(w(j,i).le.0.0) goto 12
            avres(i)=avres(i)+res(j,i)*w(j,i)
Cek  changed by ek to res*w*res*w::::>
            rms(i)=rms(i)+ (res(j,i)*res(j,i))*w(j,i)*w(j,i)
 12      continue
         davar1=davar1+rms(i)
         tres=tres+avres(i)
cek next statement
         if( (j2-nobswithw0) .le.1) goto 11
         rms(i)=sqrt( (rms(i)-avres(i)**2/(j2-nobswithw0) )
     &               /( (j2-nobswithw0) -1))
         avres(i)=avres(i)/(j2-nobswithw0)
 11   continue
      if(nitt.eq.0)then
         if(.not.single_turbo)then
            write(16,310) knobst
         endif
      endif
310   format(/,' Total number of observations is: ',i5,/)
      if(nitt.eq.0.and.isingle.ne.0)then
         if(.not.single_turbo)then
            write(16,*)'Number of observations with '//
     &                 'normalized weight 0.0 :', nobswithw0
            write(16,*)'knobs(i)   = ',knobst
            write(16,*)'nobswithw0 = ',nobswithw0
         endif
      endif
c
c     calculate data variance:
c
c     NVAReff is the actual number of unknowns to solve for !!!
c
      if( (knobst-nobswithw0).gt.nvareff)then
        davar1=(davar1-tres*tres/(knobst-nobswithw0))
     &          / ((knobst-nobswithw0)-nvareff)
        xmsqrs1=davar1*((knobst-nobswithw0)-nvareff)/(knobst-nobswithw0)
      else
         davar1=999.99
         xmsqrs1=999.99
      endif
c
      i f (.not.single_turbo) t h e n
      if(isingle.eq.0)then
         write(16,*)
         write(16,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write(16,*)'Events with  | AVRES | > 1.0 SEC are suspicious !'
         write(16,*)
         j=0
         do i=1,legs
            if(abs(avres(i)).gt.1.0)then
               write(16,'('' Event# '',i3,'' >>> '',1x,3i2.2,1x,2i2.2,
     &                    ''  AVRES ='',f6.2,'' NOBS ='',i3)')
     &                    i,iyr(i),imo(i),iday(i),ihr(i),imin(i),
     &                    avres(i),knobs(i)
               j=j+1
            endif
         enddo
         if(j.gt.0)then
         write(16,*)
         write(16,*)'^^^^^^^^^ C H E C K   these events above ^^^^^^^^'
         write(16,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         else
            write(16,*)'ZERO events of this kind found! (lucky guy!!!)'
         endif
         write(16,*)
      endif
      e n d i f
c
      return
      end ! of subr. rmsdatvar
c
      subroutine AVRESISTATIST
c
c     compute average residual statistics for the iteration just finished
c
      implicit none
      integer ifirstrun
      integer nrtotres
      real abtotres,avtotres,proz,oldres
c
      include 'vel_com.f'
c
      save oldres, ifirstrun
c
      write(16,*)
      if(nrotheres.gt.0) avotheres=avotheres/nrotheres
      if(nrotheres.gt.0) abotheres=abotheres/nrotheres
c
      if(nrrefrres.gt.0) avrefrres=avrefrres/nrrefrres
      if(nrrefrres.gt.0) abrefrres=abrefrres/nrrefrres
c
      if(nrreflres.gt.0) avreflres=avreflres/nrreflres
      if(nrreflres.gt.0) abreflres=abreflres/nrreflres
c
      nrtotres=nrotheres+nrrefrres+nrreflres
c
      if(nrtotres.gt.0)abtotres=(abotheres*nrotheres+
     &   abrefrres*nrrefrres+abreflres*nrreflres)/nrtotres
c
      if(nrtotres.gt.0)avtotres=(avotheres*nrotheres+
     &   avrefrres*nrrefrres+avreflres*nrreflres)/nrtotres
c
      write(16,'(1x,''After'',i3,'' iterations we got:'')') nitt
      write(16,*)'Average absolute & unweighted [and mean] residual of'
      write(16,'(1x,i5,'' straight and direct rays ='',f9.5,'' ['',f9.5,
     &           '']'')')nrotheres,abotheres,avotheres
      write(16,'(1x,i5,'' refracted           rays ='',f9.5,'' ['',f9.5,
     &           '']'')')nrrefrres,abrefrres,avrefrres
      write(16,'(1x,i5,'' reflected           rays ='',f9.5,'' ['',f9.5,
     &           '']'')')nrreflres,abreflres,avreflres
      write(16,*)
      if(ifirstrun.ne.10000001)then ! first run; no 'oldres' available...
         ifirstrun=10000001
         proz=0.0
      else ! NOT the first time in this routine; proz may be calculated
         if(ABS(oldres).gt.1.0e-10)then
            proz=100.*(abtotres-oldres)/oldres
         endif
      endif
      write(16,'(1x,i5,'' ALL                 RAYS ='',f9.5,'' ['',f9.5,
     &           '']'',5x,f7.2,'' %'')')nrtotres,abtotres,avtotres,proz
      oldres=abtotres
      write(16,*)
c
      return
      end ! of subr. avresistatist
c
      subroutine GAPCALC(i)
c
c     determine GAP for one event
c
      implicit none
      integer i
      include 'vel_com.f'
      integer nofgaps,j,ig
      real xstn,ystn,xhyp,yhyp,dxstnhyp,dystnhyp
c
      integer iga(200)     !  max. 200 obs. for a single event
c
c     event: i  knobs(i)  +/-e(2,i)=x    e(3,i)=y    x=+/-x(nsta,1)  y=...,2)
c     station# are stored in array  ISTM(iobs,ievent)
c
c---- compute GAP for this event
      if(knobs(i).gt.200)then
         i f (.not.single_turbo)then
         write(16,*)' WARNING:'
         write(16,*)' Event# ',i,'   Nobs = ',knobs(i),' > 200 !!!'
         write(16,*)' Array IGA(100) is too small; redimension it in'
         write(16,*)' subr. GAPCALC'
         e n d i f
         write(6,*)' Event# ',i,'   Nobs = ',knobs(i),' > 200 !!!'
         write(6,*)' Array IGA(200) is too small; redimension it in'
         write(6,*)' subr. GAPCALC'
         stop'subr. GAPCALC >>> array IGA is too small !'
      endif
      nofgaps=0
      do j=1,knobs(i)
         if(w(j,i).gt.0.0)then
            nofgaps=nofgaps+1
            xstn=x(istm(j,i),1)
            ystn=x(istm(j,i),2)
            xhyp=e(2,i)
            yhyp=e(3,i)
cek    avoiding call to atan2(zero,zero)
cek              write(6,*) ' Gapcalc: i,nofgaps,xstn,ystn,xhyp,yhyp'
cek              write(6,'(1x,i4,2x,i3,2x,4f10.3)') i,nofgaps,
cek     +                  xstn,ystn,xhyp,yhyp
            dxstnhyp=abs(xstn-xhyp)
            dystnhyp=abs(ystn-yhyp)
            if(dxstnhyp.gt.0.0001.or.dystnhyp.gt.0.0001) then
              iga(nofgaps)=57.296*ATAN2(xstn-xhyp,ystn-yhyp)
            else
              iga(nofgaps)=359
            endif
cek
            if(iga(nofgaps).lt.0) iga(nofgaps)=iga(nofgaps)+360
         endif
      enddo
      if(nofgaps.gt.0)then
         call SORTI(iga,nofgaps)
      else
         write(6,*)'WARNING: Event-# :',i,'has zero observations!'
         if(.not.single_turbo)then
            write(16,*)'WARNING: Event-# :',i,'has zero observations!'
         endif
      endif
      igap(i)=0
      ig=iga(1)-iga(nofgaps)
      if(ig.lt.0) ig=ig+360
      if(ig.gt.igap(i)) igap(i)=ig
      do j=2,nofgaps
         ig=iga(j)-iga(j-1)
         if(ig.lt.0) ig=ig+360
         if(ig.gt.igap(i)) igap(i)=ig
      enddo
c---- IGAP(i) is the gap of this event (nr. i)
c
      return
      end ! of subr. gapcalc
c
      subroutine BACKUP
c
c     go back in direction of the previous solution by doing
c     only half of the adjustments calculated in the last iteration
c
      implicit none
      include 'vel_com.f'
      integer jjj,i,n,k,iccc,k2,nl,j1,j2,ifl,k1,ksta1,j,kk1,m
      integer ksta2
      real zzz,avelo
c
      real cc(ist)
      character*1 reflch
c
      jjj=0
      do 2 i=1,legs
c---------- hypocenter backup
      n=4
c
      if(i.gt.neqs) n=1
      do 2 k=1,n
      jjj=jjj+1
      b(jjj)=b(jjj)/2.0
      if(k.eq.4)then                      ! backup depth
c
c     concept of effdeltaz no longer in use!
c
c         if(effdeltaz(i).eq.0.0)then
            e(k,i)=e(k,i)-b(jjj)          ! depth not constrained
c         else
c            effdeltaz(i)=effdeltaz(i)/2.0
c            e(k,i)=e(k,i)-effdeltaz(i)           !  depth was constrained:
c         endif                                   ! go back only half of the
         if(itopo.gt.0.and.e(k,i).lt.0.0)then
c           depth above zero... is it also above surface?
            call CHTOP(-e(2,i),e(3,i),zzz,
     &                 topo1file,topo2file)
            if(e(k,i).lt.zzz) e(k,i)=zzz  ! set hypocenter down to surface !
            if(ifixsolution.gt.0) e(k,i)=zzz ! fix depth to min_depth allowed!!
         endif
      else                                       ! effective depth-change !
            e(k,i)=e(k,i)-b(jjj)
      endif
   2  continue
ccc      iccc=mod((nitt-1),invertratio)
      iccc=mod(nitt,invertratio)
      if(iccc.ne.0) goto 510
c---- velocity readjustments for each model
      do 26 k2=1,nmod
      nl=nplay(k2)
      j1=4*neqs+nshot+laysum(k2)
      j2=j1-1+nl
      k=0
      if(.not.single_turbo)then
         write(16,3)
3        format(1x,'Velocity readjustments:')
         write(16,27) k2
27       format(1x,'Velocity model',i4)
      endif
      do 4 jjj=j1,j2
      k=k+1
      b(jjj)=b(jjj)/2.
      vp(k2,k)=vp(k2,k)-b(jjj)
      if(lowveloclay.eq.0)then
         if(k.gt.1)then
            if(vp(k2,k).lt.vp(k2,k-1))then
               i f (.not.single_turbo)then
               write(16,*)'WARNING:'
               write(16,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') k
               write(16,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               e n d i f
               write(6,'(1x,''WARNING: Tried to introduce a '',
     &                    ''low-velocity-layer! (Layer '',i2,'')'')') k
               write(6,'(1x,''Setting DVP from '',f5.2,'' to 0.0 and'',
     &                    '' VP to vp(layer_above)+0.001'')') b(jjj)
               write(6,*)
               b(jjj)=0.0
               vp(k2,k)=vp(k2,k-1)+0.001
            endif
         endif
      endif
5     format(1x,3f7.3,3x,a1)
      if(k.eq.ireflector)reflch=reflchar
      if(k.ne.ireflector)reflch=' '
      if(.not.single_turbo)then
         write(16,5) vp(k2,k),b(jjj),hp(k2,k),reflch
      endif
   4  continue
c
c    calculate and print average velocities of the model k2 :
c
      i f (.not.single_turbo) t h e n
      ifl=1
      write(16,*)
      write(16,*)'Calculation of average velocity starts at layer # ',
     &           ifl
      avelo=0
      do k=ifl+1,nplay(k2)
         avelo=avelo + ( hp(k2,k)-hp(k2,k-1) ) * vp(k2,k-1)
         write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &             ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &       hp(k2,k-1),hp(k2,k),vp(k2,k-1),avelo/(hp(k2,k)-hp(k2,ifl)),
     &       hp(k2,k)
      enddo
      write(16,*)
c
      ifl=2
      write(16,*)
      write(16,*)'Calculation of average velocity starts at layer # ',
     &           ifl
      avelo=0
      do k=ifl+1,nplay(k2)
         avelo=avelo + ( hp(k2,k)-hp(k2,k-1) ) * vp(k2,k-1)
         write(16,'('' z ='',f6.2,'' ...'',f6.2,''    v = '',f4.2,
     &             ''    AVEL = '',f4.2,'' downto z = '',f6.2)')
     &       hp(k2,k-1),hp(k2,k),vp(k2,k-1),avelo/(hp(k2,k)-hp(k2,ifl)),
     &       hp(k2,k)
      enddo
      write(16,*)
      write(16,*)
      e n d i f
c
26    continue
510   if(nsinv.eq.0) goto 511
ccc      if(mod((nitt-1),invertratio).ne.0) goto 511
      if(mod(nitt,invertratio).ne.0) goto 511
      k1=4*neqs+nshot+nltot+1
      ksta1=ksta
      if(nsp.eq.2) ksta1=(ksta/2)
      k2=k1+ksta1-1
      do 60 j=k1,k2
60    b(j)=b(j)/2.
      do 7 j=1,nsta
      cc(j)=0.0
      if(map1(j).eq.0) goto 7
      if(map1(j).gt.ksta1) goto 7
      kk1=k1-1+map1(j)
      cc(j)=b(kk1)
      ptcor(j)=ptcor(j)-cc(j)
7     continue
      if(.not.single_turbo)then
         write(16,8) (stn(m),ptcor(m),cc(m),m=1,nsta)
8        format(5(2x,a6,2f7.3))
 6       format(1x,'P correction readjustments:')
         write(16,*)
         write(16,*)'Half adjustments made'
      endif
      if(nsp.ne.2) goto 29
      k1=4*neqs+nshot+nltot+ksta1+1
      ksta2=ksta-ksta1
      k2=k1+ksta2-1
      do 61 j=k1,k2
61    b(j)=b(j)/2.
      if(.not.single_turbo)then
         write(16,39)
39       format(1x,'S correction readjustments:')
      endif
      do 28 j=1,nsta
      cc(j)=0.0
      if(map1(j).eq.0) goto 28
      if(map1(j).gt.ksta2) goto 28
      kk1=k1-1+map1(j)
      cc(j)=b(kk1)
      stcor(j)=stcor(j)-cc(j)
28    continue
      if(.not.single_turbo)then
         write(16,8) (stn(m),stcor(m),cc(m),m=1,nsta)
      endif
511   continue
      if(.not.single_turbo)then
         write(16,*)
         write(16,*)'Half adjustments made'
      endif
29    continue
c
      return
      end ! of subr. backup
c
      subroutine STEPLENGTHDAMP(damp)
c
c     set step length damping to model-vector in case of LARGE VELOCITY CHANGE
c     (  | deltavelocity |   >   MAXveladjALLOWED  )
c
      implicit none
      real damp
      include 'vel_com.f'
      integer i,j11,j22,jjj
      real btemp
c
      damp=1.0
      if(icount.eq.1) goto 900
      if(scale(6).eq.0.0) goto 900
      do i=1,nmod
         j11=4*neqs+nshot+laysum(i)
         j22=j11+nplay(i)-1
         do 121 jjj=j11,j22
            if(veladj.gt.abs(b(jjj))) goto 121
            btemp=veladj/abs(b(jjj)) ! ABS[b(jjj)] > veladj --> damping enlarged
            if(btemp.lt.damp) damp=btemp ! eff. damping=smallest of all dampings
 121     continue
      enddo
c
c     apply step length damping just calculated to ALL unknowns !!! :
c
      do jjj=1,nvar
         b(jjj)=b(jjj)*damp
      enddo
900   continue
c
      return
      end ! of subr. steplengthdamp
c
      subroutine STEPLENGTHCALC
c
c     compute the step-length STEPLEN actually applied;
c     STEPLEN is the euclidean length of the model-vector
c
      implicit none
      include 'vel_com.f'
      integer i
c
cu      real csum(4,2)
c
      steplen=0
cu      bsum=0
cu      k=0
cu      do j=1,4
cu         csum(j,1)=0.0
cu      enddo
      do i=1,nvar
cu        k=k+1
cu        if(i.le.4*neqs) csum(k,1)=csum(k,1)+b(i) ! CSUM is never used...!!??!!
cu        if(nitt.eq.1) csum(k,2)=csum(k,1)
cu        if(k.eq.4) k=0
cu        if(i.gt.4*neqs+nshot) bsum=bsum+b(i)  ! BSUM is never used... !!??!!
         steplen=steplen+b(i)*b(i)
      enddo
      steplen=sqrt(steplen)
      if(ibackups.gt.0) steplen=-steplen
      if(.not.single_turbo)then
         write(16,*)
         write(16,53) steplen
 53      format (' (Applied) Step length = ', f7.3/)
      endif
c
      return
      end ! of subr. steplengthcalc
c
      subroutine RESOLCOVAR(davari)
c
c     calculate resolution and covariance matrices.
c
      implicit none
      real davari
      include 'vel_com.f'
c
      integer mef,nef,n1,n2,n3,k,i,l,j
      real scale1,spread1,spread2,size,avresol
      real Rdiag(inva)
c
      character*80 pcard
c
c---- resolution and covariance calculations
      if(.not.single_turbo)then
         write(16,'(/////)')
         write(16,40)
   40    format(' Resolution and covariance calculations:')
         if(isingle.ne.0)then
cc       write(16,*)'CHOLESKY-decomposition:'
         write(16,*)'    RESOLUTION-matrix                   ',
     &              '    COVARIANCE-matrix'
         write(16,*)
         endif
      endif
   25 mef=4*neqs+nshot-4
      nef=4+nltot
      if(isingle.ne.0) NEF=4     ! only first 4 rows and first 4 standard-devs.
      n1=4*neqs
      n2=n1+nshot
      n3=n2+nltot
      if(iresolcalc.eq.2.and.isingle.eq.1) 
     &write(80,*)'Covariance matrice'
      do 15 k=1,nef
      do 16 i=1,nvar
   16 rht(i)=gcopy(k,i)
      l=k
      if(k.gt.4) l=k+mef
      call LUELMP(g,rht,nvar,rht)
      if(isingle.eq.0)then
         if(.not.single_turbo)then
            write(16,17) l
 17         format(1x,'resolution row',i3)
            write(16,18)  (rht(j),j=1,nvar)
         endif
         Rdiag(k)=rht(l)
      else
         write(pcard,18) (rht(j),j=1,nvar)
         do j=1,4
            Rc(l,j)=rht(j)     !  resolution-matrix
         enddo
      endif
 18   format(8f10.4)
      call LUELMP(g,rht,nvar,rht)
c---- put covariance into proper unit
      j=mod(l,4)
      if(j.eq.0) j=4
      scale1=scale(j)
      if(l.gt.n1) scale1=scale(1)
      if(l.gt.n2) scale1=scale(6)
      if(l.gt.n3) scale1=scale(5)
      do 19 j=1,nvar
ccc   19 rht(j)=rht(j)*davari*scale1
   19 rht(j)=rht(j)*scale1  ! <== UNIT covariance matrix !!!
       call FIXUNT(rht,neqs,nshot,nltot,ksta,scale,
     &             vdamp,itotmodels,inltot,nplay(1))
      if(isingle.eq.0)then
         if(.not.single_turbo)then
            write(16,20) l
 20         format(1x,'covariance row',i3)
            write(16,18)  (rht(j),j=1,nvar)
         endif
      else
         write(pcard(41:80),18) (rht(j),j=1,nvar)
         if(.not.single_turbo) write(16,*) pcard
         do j=1,4
            COVc(l,j)=rht(j)
         enddo
      endif
c      s(k)=sqrt(abs(rht(l)))  ! compiler produces warning with range checking
      s(k)=abs(rht(l))
      s(k)=sqrt(s(k))
 15   continue
c
      if(isingle.ne.0)then
         if(.not.single_turbo) write(16,*)
         call SPREADd(Rc,4,spread1)
         spread=spread1
         if(ifixsolution.eq.1)then
            do j=1,3
               do k=1,3
                  R3(j,k)=Rc(j,k)
               enddo
            enddo
            call SPREADd(R3,3,spread)
         endif
         if(ifixsolution.eq.9)then
            call SPREADd(Rc(1,1),1,spread)
         endif
         call SPREADb(Rc,4,spread2)
         size=0.0
         do j=1,4
            size=size+COVc(j,j)
         enddo
         if(.not.single_turbo)then
            write(16,'(1x,''D-Spread(R) = '',f6.3,
     &                    ''   B-G Spread(R) = '',f6.3,
     &                 ''       Size (C) = '',f6.3
     &             )') spread,spread2, size
         endif
      endif
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,21)
 21      format(' standard deviation of selected model parameters:')
      endif
      if(isingle.ne.0)then
        avresol=0.0
        do j=1,nef
           avresol=avresol+Rc(j,j)
        enddo
        avresol=avresol/float(nef)
        i f (.not.single_turbo) t h e n
        write(16,*)
        write(16,*)'   OT (sec)   X (km)    Y (km)    Z (km) '
        e n d i f
        write(6,'(23x,''  OT (sec)   X (km)    Y (km)    Z (km) '')')
c        write(6,'(1x,''CHOLESKY  (othet='',f5.3,'') :'')') othet
        write(6,'('' Sigma (CHD):         '',4f10.4)') (s(j),j=1,nef)
        write(6,'('' Resolution (CHD):    '',4f10.4,x,
     &            ''D-spread ='',f6.3)') (Rc(j,j),j=1,nef) , spread
        write(6,'('' Data Variance      = '',f10.4)') davari
ccc        call SPREADd(Rs,4,spread3)
ccc        do j=1,4
ccc           COVs(j)=SQRT(COVs(j))  !  <-- standard deviation
ccc        enddo
        write(6,'('' Singular values:     '',4f10.4,5x,''ALE ='',f7.3)')
     &               (SV(j),j=1,nef) , ale(1)
ccc        write(6,'('' Sigma (SVD):         '',4f10.4)') (COVs(j),j=1,nef)
ccc        write(6,'('' Resolution (SVD):    '',4f10.4,3x,
ccc     &            ''D-spread ='',f6.3)') (Rs(j,j),j=1,nef) , spread3
      endif
      if(.not.single_turbo) write(16,22) (s(j),j=1,nef)
 22   format(8f10.4)
      if(.not.single_turbo)then
         write(16,*)
         write(16,*)'Rdiag of selected model parameters:'
         write(16,22) (Rdiag(j),j=1,nef)
         write(16,*)
         write(16,*)
      endif
c
      return
      end ! of subr. resolcovar
c
      subroutine TRAVDERIV(raytype,nl,mll,v1,vsq1,
     &                     rp,nrp,x2,y2,z2,ss,r1,r2,ievent,inobs)
c
c     compute traveltime derivatives with respect to all the unknowns inverted
c     for !
c
      implicit none
      include 'vel_com.f'
c
      character*(*) raytype
      integer nl,mll,nrp,ievent,inobs,ii,j,jx,jndex,l,jb
      integer idownward,i
      real rp(3,inrpmax), f(inltot),zmax,r1,r2,dtdd
      real v1(nl),vsq1(nl),x2(nrp),y2(nrp),z2(nrp),ss(nrp)
c
      if(isingle.ne.0.and.iturbo.eq.1)then
         continue
      else
         do 2 ii=1,nl
2        dtdv(ii)=0.
      endif
      do 3 ii=1,3
3     dtdr(ii)=0.
c
      if(ifixsolution.eq.9) RETURN   ! hold LAT/LON/DEPTH fixed !!!
c
c     No trvdrv's for reseted obs:
c
      if(isingle.ne.0.and.w(inobs,ievent).eq.0.0) RETURN
c
      if(raytype.eq.'direct')    goto 5
      if(raytype.eq.'refracted') goto 6
      if(raytype.eq.'reflected') goto 8
      stop'TRAVDERIV>>> illegal raytype!!!'
c
c   direct ray
c
5     continue
      do 9 j=1,nrp-1
      jx=nrp-j
      x2(j)=rp(1,j+1)-rp(1,j)
      y2(j)=rp(2,j+1)-rp(2,j)
      z2(j)=rp(3,j+1)-rp(3,j)
      ss(j)=sqrt(x2(j)**2 + y2(j)**2 + z2(j)**2)
      if(isingle.ne.0.and.iturbo.eq.1)then
         continue
      else
         dtdv(jx)= -ss(j)/vsq1(jx)
      endif
9     continue
      dtdr(1)=-x2(1)/(v(jl)*ss(1))
      dtdr(2)=-y2(1)/(v(jl)*ss(1))
      dtdr(3)=-z2(1)/(v(jl)*ss(1))
      goto 40
c
c    refracted first path
c
6     continue
      do 171 j=1,nrp
      z2(j)=rp(3,j)
171   continue
      call MAXRI(nrp,z2,zmax,jndex) ! determine MAX (=zmax) of z2 (vertical
c                                     path-lengths in each layer)
      l=nl
12    if(h(l).le.(zmax+.01)) goto 11
      l=l-1
      goto 12
11    continue
      jb=l
c
      dtdd=1/v(jb)
      dtdr(1)=(r1/delta)*dtdd
      dtdr(2)=(r2/delta)*dtdd
      dtdr(3)=-sqrt(vsq(jb)-vsq(jl))/(v(jb)*v(jl))
c
      if(isingle.ne.0.and.iturbo.eq.1) goto 40  ! do not calc. dt/dv  !!!
c
c     thk(j) -- thickness of layer j
c     jl -- event layer
c     jb -- bottomming layer
c     tkj -- depth of event within event layer
c
      do 201 j=1,nl
      f(j)=1.
      if(j.ge.jl) f(j)=2.
201   if(j.gt.jb) f(j)=0.
      do 202 j=1,jb-1
      dtdv(jb)=dtdv(jb)+ thk(j)*v1(j)*f(j)/(vsq1(jb)*
     &sqrt(vsq1(jb)-vsq1(j)))
      dtdv(j)=-thk(j)*v1(jb)*f(j)/(vsq1(j)*
     &sqrt(vsq1(jb)-vsq1(j)))
202   continue
      dtdv(jl)=dtdv(jl)+tkj*v1(jb)/(vsq1(jl)*
     &sqrt(vsq1(jb)-vsq1(jl)))
      dtdv(jb)=dtdv(jb)-tkj*v1(jl)/(vsq1(jb)*
     &sqrt(vsq1(jb)-vsq1(jl))) -delta/vsq1(jb)
c
      goto 40
c
c     reflected phase
c
  8   continue
c  avoid hypocenter on layer-boundary:
chrm      if(rp(1,1).eq.rp(1,2).and.
chrm     &   rp(2,1).eq.rp(2,2).and.
chrm     &   rp(3,1).eq.rp(3,2))then
chrm            rp(3,1)=rp(3,1)-0.0001   ! move hypocenter 10 cm up
chrm      endif
c
ctest      dtdr(1)=DTDDrefl*(r1/delta)      ! valid for REFLECT
ctest      dtdr(2)=DTDDrefl*(r2/delta)
ctest      dtdr(3)=DTDHrefl
c
      jx=0
      do j=1,nrp-1
         x2(j)=rp(1,j+1)-rp(1,j)
         y2(j)=rp(2,j+1)-rp(2,j)
         z2(j)=rp(3,j+1)-rp(3,j)
         ss(j)=sqrt(x2(j)**2 + y2(j)**2 + z2(j)**2)
         if(j.eq.1)then
            dtdr(1)=-x2(1)/(v(jl)*ss(1))          ! valid for REFLECT1
            dtdr(2)=-y2(1)/(v(jl)*ss(1))
            dtdr(3)=-z2(1)/(v(jl)*ss(1))
         endif
         if(isingle.ne.0.and.iturbo.eq.1) goto 40 ! do not calc. dt/dv  !!!
         idownward=0
         if(z2(j).gt.0.0) idownward=1
         if(j.eq.(nrp-1).and.idownward.eq.1) idownward=0
         if(idownward.eq.1)then  !  DOWNwards
            jx=jl+j-1
c            if(vsq1(jx).le.0)then
c               write(6,*)'jx=',jx,'jl=',jl,'tkj=',tkj
c               write(6,*)'j=',j,'vsq1(jx)=',vsq1(jx),'nrp=',nrp
c            endif
            dtdv(jx)=-ss(j)/vsq1(jx)
         else   !  UPward
            if(jx.eq.mll)then
               dtdv(mll)=dtdv(mll)-ss(j)/vsq1(mll)
               jx=jx-1
            else
               dtdv(jx)=dtdv(jx)-ss(j)/vsq1(jx)
               jx=jx-1
            endif
         endif
      enddo
      if(jx.ne.0)stop'TRAVDERIV>>> did not reach the top...!'
      goto 40
c
c
40    continue
c
      if(ifixsolution.eq.1) dtdr(3)=0.0   ! hold depth fixed !!!
c
cek   test if nrp gt. inrpmax:::>
      if(nrp.gt.inrpmax) then
       if(.not.single_turbo)then
          write(16,2000) i,inobs,nrp
 2000  format(//,2x,'travderiv>>>, label 2000:  ',
     & 'nrp greater than inrpmax:',/,
     & 5x,'event nr=',i6,2x,'nr. of obs=',i4,5x,4hnrp=,i4,/,
     & 3x,'program stops here',/)
       endif
       stop'subr. TRAVDERIV >>> nrp > nrpmax'
      endif
c
      return
      end ! of subr. travderiv
c
      subroutine REJECTOBS(i,nobs,iresflag)
c
c     reject the current observation for further calculation in case the
c     reading (or the seismologist who read it...) seems to be nuts:
c
      implicit none
      integer i,nobs,iresflag
      include 'vel_com.f'
c
      real wsum,xkndw
      integer knobst,iii,iobswt
      if(w(nobs,i).eq.0.0) RETURN
      if( (knobs(i)-nobswithw0) .eq. nvar ) RETURN
c
cc         w(nobs,i)=0.0039  !  = 1./(2**(2*reading_weight)) = 1./256
c
c     re-normalize weights of this event :
c
      nobswithw0=nobswithw0+1
      wsum=0.0
      knobst=knobs(i)
      do iii=1,knobst
         iobswt=kpwt(iii,i)
         if(iobswt.lt.4.and.w(iii,i).ne.0.0)then
            if(sphase(iii,i).eq.1.0.or.sphase(iii,i).eq.2.0)then
               w(iii,i)=swtfac*1.0/(2**(iobswt*2))   ! S-phase or s-p phase
            else
               w(iii,i)=       1.0/(2**(iobswt*2))   ! P- or M-phase
            endif
         else
            w(iii,i)=0.0   ! observation-weight 4 ==> don't use this arrival
         endif
         if(iii.eq.nobs) w(iii,i)=0.0
         wsum=wsum+w(iii,i)
      enddo
      xkndw=float(knobst-nobswithw0)/wsum
      do iii=1,knobst
         w(iii,i)=w(iii,i)*xkndw
      enddo
      if(.not.single_turbo)then
      write(16,*)'WARNING:'
      write(16,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &   '' ='',f7.2,''; ABS > 2.0 ---> weight set to zero !'')')
     &   nitt,isingle,nobs,res(nobs,i)
      endif
      write(6,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &   '' ='',f7.2,''; ABS > 2.0 ---> weight set to zero !'')')
     &   nitt,isingle,nobs,res(nobs,i)
      if(.not.single_turbo)then
         write(16,*)'knobs(i)   = ',knobst
         write(16,*)'nobswithw0 = ',nobswithw0
      endif
cc      write(6,*)'knobs(i)   = ',knobst
cc      write(6,*)'nobswithw0 = ',nobswithw0
cc      write(6,*)
      iresflag=1  ! inhibit a pre-stopping in subr. OUTPUT due to
                  ! increased datvar; give location a second chance.
c
      return
      end ! of subr. rejectobs
c
      subroutine REVIVEOBS(i,nobs,iresflag)
c
c     If observationweight < 4 and weight=0.0 then set 'original' weight
c     if residual has become smaller since it has been reseted:
c
      implicit none
      integer i,nobs,iresflag
      include 'vel_com.f'
c
      real wsum,xkndw
      integer knobst,iii,iobswt
      if(w(nobs,i).ne.0.0) RETURN
      if(kpwt(nobs,i).eq.4) RETURN
c
cc         w(nobs,i)=0.0039  !  = 1./(2**(2*reading_weight)) = 1./256
c
c     re-normalize weights of this event :
c
      nobswithw0=nobswithw0-1
      wsum=0.0
      knobst=knobs(i)
      do iii=1,knobst
         iobswt=kpwt(iii,i)
         if(iii.ne.nobs.and.w(iii,i).eq.0.0)then
            w(iii,i)=0.0   ! keep zero-weights as they are !!!
         else
            if(iobswt.lt.4)then
               if(sphase(iii,i).eq.1.0.or.sphase(iii,i).eq.2.0)then
                  w(iii,i)=swtfac*1.0/(2**(iobswt*2))   ! S-phase
               else
                  w(iii,i)=       1.0/(2**(iobswt*2))   ! P- or M-phase
               endif
            else
               w(iii,i)=0.0   ! observation-weight 4 ==> don't use this arrival
            endif
         endif
         wsum=wsum+w(iii,i)
      enddo
      xkndw=float(knobst-nobswithw0)/wsum
      do iii=1,knobst
         w(iii,i)=w(iii,i)*xkndw
      enddo
      if(.not.single_turbo)then
      write(16,*)'WARNING:'
      write(16,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &           '' ='',f7.2,''; ABS < 1.0 ---> weight revived !'')')
     &           nitt,isingle,nobs,res(nobs,i)
      endif
      write(6,'('' Iteration# '',i3,'' Event#'',i4,'' Res#'',i3,
     &           '' ='',f7.2,''; ABS < 1.0 ---> weight revived !'')')
     &           nitt,isingle,nobs,res(nobs,i)
      if(.not.single_turbo)then
         write(16,*)'knobs(i)   = ',knobst
         write(16,*)'nobswithw0 = ',nobswithw0
      endif
cc      write(6,*)'knobs(i)   = ',knobst
cc      write(6,*)'nobswithw0 = ',nobswithw0
cc      write(6,*)
      iresflag=1  ! inhibit a pre-stopping in subr. OUTPUT due to
                  ! increased datvar; give location a second chance.
c
      return
      end ! of subr. reviveobs
c
      subroutine DETNOFUNKNOWNS
c
      implicit none
      include 'vel_com.f'
      integer lipeff
c
c     determine the number of unknowns NVAR and
c     calculate the number of equations LIP to be solved :
c
      if(icount.ne.0)then
         lip=4*neqs + nshot +1 ! do NOT invert for velocity parameters
         lipeff=lip
      else
         lip=4*neqs + nshot + nltot + ksta +1 ! invert for velocity parameters
         lipeff=4*neqs + nshot + nltot + nstaeff +1 ! invert for vel. parameters
         if(nsinv.eq.0) lip=4*neqs + nshot + nltot + 1 !do NOT inv. for sta-corr
      endif
c     NVAR is the number of unknowns to solve for :
      nvar=lip-1
c     NVAREFF is the REAL number of unknowns to solve for!!!
      nvareff=lipeff-1
c     KVAR is the number of elements on or above the main diagonal of G=(At(A))
      kvar=nvar*lip/2              !  cuk !BELOW!! cuk
c
c
      if(.not.single_turbo)then
         write(16,*)
         write(16,'('' Number of unknowns (for array-indexing): '',
     &              '' nvar = '',i4)') nvar
         write(16,'('' Number of effective unknowns        : '',
     &              '' nvareff = '',i4)') nvareff
         write(16,'('' Number of elements on/below main diagonal '',
     &              ''of matrix G = At*A : kvar = '',i7)') kvar
         write(16,*)
      endif
c
      return
      end ! of subr. detnofunknowns
c
      subroutine ACTUALSTATIONS
c
cek added separate report of P and S phases 20.8.98
c
      implicit none
      include 'vel_com.f'
      integer i,nofreadings,k,nobsp(ist),nobss(ist),nofreadp,nofreads
c
      do i=1,nsta
         nactualsta(i)=0
         nobsp(i)=0
         nobss(i)=0
      enddo
      nstaeff=0
      nofreadings=0
      do i=1,legs  ! = neqs+nshot
         do k=1,knobs(i)
            nactualsta( istm(k,i) ) = nactualsta( istm(k,i) ) + 1
            if(nsp.eq.2) then
             if(sphase(k,i).eq.0.0) nobsp(istm(k,i))=nobsp(istm(k,i))+1
             if(sphase(k,i).eq.1.0) nobss(istm(k,i))=nobss(istm(k,i))+1
            endif
            nofreadings=nofreadings+1
         enddo
      enddo
      if(.not.single_turbo) write(16,*)
      do i=1,nsta
         if(nactualsta(i).gt.0)then
            if(.not.single_turbo)then
               if(nsp.eq.2) then
               write(16,'('' readings for station '',a6,'' : tot='',
     &                  i4,''  P:'',i4,''  S:'',i4)')
     &                  stn(i),nactualsta(i),nobsp(i),nobss(i)
               else
               write(16,'('' readings for station '',a6,'' :'',i4)')
     &                  stn(i),nactualsta(i)
               endif
            endif
            nstaeff=nstaeff+1
         endif
         nofreadp=nofreadp+nobsp(i)
         nofreads=nofreads+nobss(i)
      enddo
      if(.not.single_turbo)then
         write(16,*)
         write(16,'('' Total number of stations with readings:'',i4)')
     &           nstaeff
         write(16,*)
         write(16,'('' Total number of readings: '',i7)') nofreadings
         write(16,'('' Total number of P readings: '',i7)') nofreadp
         write(16,'('' Total number of S readings: '',i7)') nofreads
         write(16,*)
      endif
c
      RETURN
      end ! of subr. actualstations
c
c
      subroutine FINALSTARESI
c
c     computes final station residuals for all stations used in this run.
c
      implicit none
      include 'vel_com.f'
c
      integer i,j,m,k,iz,iwarn
      real stcor1
      real del(ist)
chrm      real phz(3)
      character*1 phz(3)
      real aa(ist),bb(ist),dd(ist),ee(ist)
      integer icc(ist), iccs(ist)
      real aas(ist),bbs(ist),dds(ist),ees(ist)
      character*1 cns,cew,reflch
      character*40 titl
      character cline*80, sta*6
      character*1 phzz(ist)
      data phz(1),phz(2),phz(3)/'S','P','m'/
c
      if(isingle.eq.0)then
         if(iturbo.eq.0)then
            write(16,*)
            write(16,*)
            write(16,*)
            write(16,*)
         else
            write(16,*)'TURBO-option is set; residuals are NOT printed '
     &                 //'for each event!'
         endif
      endif
      write(16,*)
      if(isingle.eq.0)then
         if(iturbo.eq.0)then
            write(16,1)
 1          format(1h1)
         else
            goto 8001    ! do NOT output residuals for each event if TURBO set.
         endif
      endif
c
c     output residuals for each event seperately:
c
      do 2 i=1,legs
         write(16,*)
         if(imin(i).lt.0)then  !  U.K.  3.Feb.87
            imin(i)=imin(i)+60
            ihr(i)=ihr(i)-1
         endif
         if(isingle.eq.0)then
          write(16,'(1x,''Station residuals for event='',i4,
     &               3x,3i2.2,1x,2i2.2,1x,f5.2)')
     &          i,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i)
         else
          write(16,'(1x,''Station residuals for event='',i4,
     &               3x,3i2.2,1x,2i2.2,1x,f5.2)')
     &          isingle,iyr(i),imo(i),iday(i),ihr(i),imin(i),e(1,i)
         endif
c        no more active
c3       format(1h ,'Station residuals for event=',i4,
c    &               3x,3i2.2,1x,2i2.2,1x,f5.2)
         write(16,51)
51       format(1x,2('sta ph wt  res   ttime delta',5x))
         do 52 j=1,knobs(i)
            del(j)=sqrt((e(2,i)-d(j,1,i))**2 + (e(3,i)-d(j,2,i))**2)
            phzz(j)='p'
            if(sphase(j,i).eq.1.0) phzz(j)='s'
            if(sphase(j,i).eq.-1.0) phzz(j)='m'
52       continue
         write(16,53) (smn(j,i),phzz(j),kpwt(j,i),res(j,i),tctime(j,i),
     &                 del(j),j=1,knobs(i))
53       format(2(1x,a6,1x,a1,1x,i2,f7.3,2f6.2,4x))
2     continue
      write(16,*)
      if(isingle.ne.0) goto 8001
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)
 8001 continue
c
c     output station statistics (nobs, avres, ...) :
c
      if(isingle.eq.0) then
        write(16,7777) nsp
7     format(/,1x,'sta phase nobs avres  avwres    std    wsum    ',
     &         'delay',/)
        write(16,7)
7777  format(//,2x,' station statistics, remember nsp was set to:',i2,/)
      endif
48    do m=1,nsta
         aa(m)=0.0
         bb(m)=0.0
         icc(m)=0
         dd(m)=0.0
         ee(m)=0.0
         aas(m)=0.0
         bbs(m)=0.0
         iccs(m)=0
         dds(m)=0.0
         ees(m)=0.0
      enddo
c loop 9: collect residual data from all events
      do 9 m=1,nsta
         sta=stn(m)
         do 10 i=1,legs
            k=knobs(i)
            do 11 j=1,k
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.0.) goto 12    ! P
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.-1.0) goto 12  ! refl. P
11          continue
            goto 110
12          continue
            aa(m)=aa(m)+res(j,i)*w(j,i)
cek         res: *w(j,i)
            bb(m)=bb(m)+res(j,i)*res(j,i)*w(j,i)*w(j,i)
            ee(m)=ee(m)+res(j,i)
            dd(m)=dd(m)+w(j,i)
            icc(m)=icc(m)+1
110         continue
            if(nsp.eq.1) goto 10
c
chrm        Doing S and S-P phases
c
            do 101 j=1,k
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.1.) goto 102 ! s-phase
               if(sta.eq.smn(j,i).and.sphase(j,i).eq.2.) goto 102 ! s-p phase
101         continue
            goto 10
102         continue
            aas(m)=aas(m)+res(j,i)*w(j,i)
c           res: *w(j,i)
            bbs(m)=bbs(m)+res(j,i)*res(j,i)*w(j,i)*w(j,i)
            ees(m)=ees(m)+res(j,i)
            dds(m)=dds(m)+w(j,i)
            iccs(m)=iccs(m)+1
10       continue
9     continue
c end loop 9
c
c     Station file output:
c
      if(istaout.gt.0)then
         open(12,file=stafile,status='unknown')
         write(12,1201) fm
 1201    format(a80)
      endif
c
c begin loop 13: print output
      do 13 m=1,nsta
c
c                              write on file 12
         if(istaout.gt.0)then
            iz=ielev(m)
            if(xla(m).lt.0.0)then
               cns='S'
               xla(m)=-xla(m)
            else
               cns='N'
            endif
            if(xlo(m).lt.0.0)then
               cew='E'
               xlo(m)=-xlo(m)
            else
               cew='W'
            endif
            write(cline,fm) stn(m),xla(m),cns,xlo(m),cew,
     &                   iz,model(m),map1(m),
     &                   ptcor(m),stcor(m)
            if(cns.eq.'S') xla(m)=-xla(m)
            if(cew.eq.'E') xlo(m)=-xlo(m)
            if(m.eq.1)then
               cline(47:)='       lon,z,model,icc,ptcor,stcor'
            endif
            write(12,'(a)') cline
         endif
c                            end writing on file 12
c
         if(dd(m).eq.0.or.icc(m).lt.2) goto 1013
         bb(m)=(bb(m)-aa(m)**2/dd(m))*icc(m)/(dd(m)*(icc(m)-1))
         iwarn=0
         if(bb(m).lt.-0.1)then
           iwarn=1
           write(16,'(5x,''WARNING: Station = '',a6,
     &          5x,''!!! Variance  bb('',i3,'') = '',f7.3,'' < 0 !!'')')
     &          stn(m),m,bb(m)
cek           write(6,'(5x,''WARNING: Station = '',a6,
cek     &         5x,''!!! Variance  bb('',i3,'') = '',f7.3,'' < 0 !!'')') 
cek     &         stn(m),m,bb(m)
cek           write(6,*)
            bb(m)=0.0
         else
            bb(m)=abs(bb(m))
         endif
         if(iwarn.eq.0) then
            bb(m)=sqrt(bb(m))
         else
         endif
         aa(m)=aa(m)/dd(m)
         ee(m)=ee(m)/icc(m)
c
cek  phz(2)=  P-phases
c
cek         if(iwarn.eq.0)then
         write(16,'(1x,a6,3x,a1,i4,3f8.4,1x,2f8.4)') 
     &             stn(m),phz(2),icc(m),ee(m),aa(m),bb(m),dd(m),ptcor(m)
cek         else
cek            write(16,'(1x,a6,3x,a1,i4,2f8.4,''   ?.?.?''1x,2f8.4,
cek     &              '' <--- !!'')')
cek     &             stn(m),phz(2),icc(m),ee(m),aa(m),dd(m),ptcor(m)
cek            iwarn=0
cek         endif
1013     if(nsp.eq.1) goto 39
cek         write(16,'(1x,a6,3x,a1,i4,25x,f8.4)') 
cek     &                         stn(m),phz(1),iccs(m),dds(m)
c  print output for S-wave data  (nsp=2 or 3):
         if(dds(m).eq.0.or.iccs(m).lt.2) goto 13
         bbs(m)=(bbs(m)-aas(m)**2/dds(m))*iccs(m)/(dds(m)*(iccs(m)-1))  
         iwarn=0
         if(bbs(m).lt.-0.1)then
           iwarn=1
           write(16,'(5x,''WARNING: Station = '',a6,
     &        5x,''!!! Variance  bbs('',i3,'') = '',f7.3,'' < 0 !!'')')
     &        stn(m),m,bbs(m)
cek           write(6,'(5x,''WARNING: Station = '',a6,
cek     &         5x,''!!! Variance  bb('',i3,'') = '',f7.3,'' < 0 !!'')') 
cek     &         stn(m),m,bb(m)
cek           write(6,*)
         else
            bbs(m)=abs(bbs(m))
         endif
         if(iwarn.eq.0) then
            bbs(m)=sqrt(bbs(m))
         else
            bbs(m)=0.0
         endif
         aas(m)=aas(m)/dds(m)
         ees(m)=ees(m)/iccs(m)
c
c   (nsp=3):
         if(nsp.eq.3) then
            stcor1=ptcor(m)*vpvs
c   (nsp=2):
         else
            stcor1=stcor(m)
         endif
c
cek  phz(1)=  S-phases
c
cek         if(iwarn.eq.0)then
         write(16,'(1x,a6,3x,a1,i4,3f8.4,1x,2f8.4)') 
     &             stn(m),phz(1),iccs(m),ees(m),aas(m),bbs(m),dds(m),
     &             stcor1
cek         else
cek            write(16,'(1x,a6,3x,a1,i4,2f8.4,''   ?.?.?''1x,2f8.4,
cek     &              '' <--- !!'')')
cek     &             stn(m),phz(1),iccs(m),ees(m),aas(m),dds(m),stcor1
cek            iwarn=0
cek         endif
c
 14      format(1x,a6,3x,a1,i4,3f8.4,1x,2f8.4)
 39      continue
 13   continue
c
c end loop print one station
c
      write(16,*)
      write(16,*)
      if(istaout.gt.0)then
         write(12,*)
         close(12)
      endif
c
chrm  Output of the final model
c
c M. Zhang changed the output format to fit the origional model format
      if(istaout.eq.2)then
         open(12,file='velout.mod',status='unknown')
         write(12,'(a40)')'Output model by isingle = 0:'
c	 write(12,*)(nplay(m),m=1,nmod)
      do m=1,nmod
            write(12,'(i3)')nplay(m)
            reflch=' '
            titl=' '
        do i=1,nplay(m)
        if(i.eq.1)then
c   The surface layer velocity is less than the lower one
            if(vp(m,i) .gt. vp(m,i+1))then
                write(12,'(f5.2,5x,f7.2,2x,f7.3,3x,a1,1x,a40)')
     &              vp(m,i+1),hp(m,i),vdamp(m,i),reflch,titl
            else
                write(12,'(f5.2,5x,f7.2,2x,f7.3,3x,a1,1x,a40)')
     &              vp(m,i),hp(m,i),vdamp(m,i),reflch,titl
            endif
        else
            write(12,'(f5.2,5x,f7.2,2x,f7.3,3x,a1)')
     &          vp(m,i),hp(m,i),vdamp(m,i),reflch
        endif
        enddo
      enddo
      close(12)
      endif 
      return
      end ! of subr. finalstaresi
c
cek   end of vel_kern.f
c
cek    begin of vel_mag.f
c
      SUBROUTINE MAGNITUDE(ievent,itime)   ! Urs Kradolfer, July 1987
c
c     calculate magnitude for an event with depth Z and NRP observations
c
      implicit none
      integer ievent,itime
      include 'vel_com.f'
c
      integer iseis,iscon,isdmp,isamp,ier,i
      real cormag,sconst,sdampf,voltgain,xmag
      character*2 ifilt
      character staname*6
c
      nmag=0
      xmagnitude=0.
      sdxmagnitude=0.
      do i=1,knobs(ievent)
         if(kpwt(i,ievent).gt.4) goto 40 ! weights > 4 are no readings!!!
c                              station# are stored in array  ISTM(iobs,ievent)
c                              smn(iobs,iev)  is the station-name
         write(staname,'(a6)') smn(i,ievent)
         call STINP(itime,staname,ifilt,iseis,iscon,
     &                            isdmp,isamp,cormag,ier)
         if(ier.lt.0) goto 40  ! seismo-parameters of station not found
         if(iyr(ievent).ge.1984.and.ifilt.eq.'DE') ifilt='AD'
         sconst=float(iscon)/10.
         sdampf= float(isdmp)/10.
         voltgain=float(isamp)
         delta=SQRT(  (e(2,ievent)-d(i,1,ievent))**2
     &              + (e(3,ievent)-d(i,2,ievent))**2  )
         xmag=0.
         call MUK(delta,e(4,ievent),ifilt,iseis,sconst,sdampf,
     &            voltgain,cormag,amx(i),prx(i),xmagni(i))
         if(xmagni(i).eq.-13.) goto 40
         nmag=nmag+1
         xmag=xmagni(i)
         xmagnitude=xmagnitude+xmag
         sdxmagnitude=sdxmagnitude+xmag**2
  40     continue
      enddo
      if(nmag.ne.0)then
         if(nmag.ge.2)then
            sdxmagnitude=
     &      SQRT( (nmag*sdxmagnitude-xmagnitude*xmagnitude)
     &           / (nmag*(nmag-1)) )
         else
             sdxmagnitude=0.
         endif
         xmagnitude=xmagnitude/nmag
      else
         xmagnitude=0.0
      endif
c
      return ! calculated magnitude is: AVXM +/- SDXM   and all xmagni(1...nobs)
c
      END ! of subr. magnitude
c
      subroutine MUK(epdist,depth,ifilter,isecpendel,seismkonst,
     &               seismdamp,voltgain,stacor,devampl,period,xmag)
c
c------------------------------- Urs Kradolfer 1984 -------------------
C
C     Berechnet die Magnitude Ml fuer Epizentraldistanz <=  700 km
C     Berechnet die Magnitude Mb fuer Epizentraldistanz > 2200 km
c
c     Input:
c             Epizentraldistanz         [km]       epdist  >=0.  {real}
c                    "                  [deg]      epdist  < 0.  {real}
c             Herdtiefe                 [km]       depth         {real}
c             filter-typ (z.b. DE, AD, PD)         filter-typ    {char}
c             Seismometer-Grenzperiode  [sec]      isecpendel    {integer}
c                  "     -Konstante     [V/(cm/s)] seismkonst    {real}
c                  "     -Daempfung     [1]        seismdamp     {real}
c             Verstaerkung (Elektronik) [1]        voltgain      {real}
c             Stationskorrektur         [Mag]      stacor        {real}
c             Max. Ablese-Amplitude     [mm]       devampl       {real}
c                 (Devco oder Plot)
c             dominierende Periode      [sec]      period        {real}
c
c     Output: Magnitude (Ml oder Mb)    [Mag]      xmag          {real}
c
c        >>   Falls die Magnitude nicht berechnet werden kann, erfolgt
c        >>   ein RETURN und der Output-Parameter xmag wird -13. gesetzt!
c
c     Zur Subr. muk gehoert ebenfalls die Subr. UFWABO und die
c     zur letzteren gehoerenden vier Complex-Functions.
c---------------------------------------------------------------------------
c
      implicit none
      real epdist,depth,seismkonst,seismdamp,voltgain
      real stacor,devampl,period,xmag
      integer isecpendel
c
      real epdistkm,ampl,waampl,boampl,delta,sigma,a,b
      integer iampltype,i
      character*2 ifilter
      real RDELT(12), SIGAR(12)
      DATA RDELT/4.0,6.,13.,16.,28.,87.,
     & 114.,120.,134.,141.,146.,180./,
     &     SIGAR/6.2,7.,7.0,5.8,6.6,7.0,7.5,
     & 7.5, 6.9, 7.1, 6.9, 6.9 /
      if(isecpendel.le.0.or.seismkonst.le.0..or.seismdamp.le.0..or.
     &   voltgain.le.0..or.devampl.le.0..or.period.le.0.)then
         xmag=-13.    ! Magnitude kann nicht berechnet werden
         return
      endif
      xmag=-13.  ! falls xmag nicht anders berechnet wird !!
      if(epdist.lt.0.)then
         epdistkm= -epdist/360.*40030     ! epdist < 0. --> [epdist] = Grad
      else
         epdistkm=sqrt(epdist**2+depth**2)
      endif
      if(epdistkm.le.2200.)iampltype=1    ! Wood-Anderson --> Ml
      if(epdistkm.gt.2200.)iampltype=2    ! Bodenbewegung --> Mb
      if(ifilter.eq.'DE'.or.ifilter.eq.'AD')
     &call ufwabo(isecpendel,seismkonst,seismdamp,voltgain,
     &            devampl,period,iampltype,ampl,ifilter)
      if(ifilter.eq.'PD')
     &call mpdr2(isecpendel,seismkonst,seismdamp,voltgain,
     &            devampl,period,iampltype,ampl)
      if(ampl.le.1.0e-10)then
         xmag=-13.  ! avoid taking LOG of non-positive value !!!
         RETURN
      endif
      if(iampltype.eq.1)waampl=ampl/2. ! p-p-Ampl. --> 0-p-Ampl.
      if(iampltype.eq.2)boampl=ampl
      if(epdistkm.gt.2200.)then
         delta=epdistkm/40030.*360.     ! delta:=epdistkm in Grad
         if(delta.gt.180.)then
            delta=180.
            sigma=sigar(12)
            goto 2
         endif
         i=0
1        i=i+1
         if(delta.gt.rdelt(i))goto 1
         ! jetzt ist delta < rdelt(i)     --> Linearisieren
         ! y=a*x+b
         a=(sigar(i)-sigar(i-1))/(rdelt(i)-rdelt(i-1))
         b=sigar(i)-a*rdelt(i)
         sigma=a*delta+b
2        xmag=log10(boampl/period) + sigma + stacor
         RETURN
      endif
      if(epdistkm.gt.0..and.epdistkm.le.60.)then
         xmag=log10(waampl) + 0.018 *epdistkm+1.77 + 0.40
      endif
      if(epdistkm.gt.60..and.epdistkm.le.700.)then
         xmag=log10(waampl) + 0.0038*epdistkm+2.62 + 0.40
      endif
      if(epdistkm.gt.1100..and.epdistkm.le.1700.)then
         xmag=log10(waampl) + 0.0029*epdistkm+3.40 + 0.40 - 2.  ! EMPIRISCH
      endif
      if(xmag.eq.-13.) RETURN     ! xmag konnte nicht berechnet werden
      xmag=xmag+stacor   ! Stationskorrektur Cs
c
      return
      end
c
      subroutine MPDR2(isecpendel,seismkonst,seismdamp,
     &               voltgain,pdrampl,period,iampltyp,ampl)
c
c------------------------------------ Kaspar G. Renggli 1984 ---------------
C
c     bodenamliptude fuer pdr2-stationen des nw-ch-netzes
c     parameter wie bei subr. muk, ausser pdrampl in [1/10 um/sec]
c---------------------------------------------------------------------------
c  
      implicit none
      integer isecpendel,iampltyp
      real seismkonst,seismdamp,voltgain,pdrampl,period,ampl
      complex G,gseis,ghp,glp2,glp2a,glp3,j,jom
      real zpi,ts,s,hs,om,om2,omg1,omg2,omg3,omg4,omg5,twa
      real omw,hw,vwa
      data zpi/6.28319/ j/(0.,1.)/
      TS=isecpendel
      S=seismkonst
      hs=seismdamp
      om=zpi/period
      jom=j*om
      om2=om**2
c-----Uebertragungsfunktionen fuer Seismometer und Filter des PDR-2-Systems---
c
c------------ Seismometer als Wegaufnehmer (nur Mechanik) -------------
      omg1=zpi/TS
      gseis=jom*om2/(omg1**2+2*jom*omg1*hs-om2)
c------------- 2 Hochpaesse 1. Ordnung  =>  2. Ordnung, fc=0.3 Hz -----
      omg2=zpi*0.3
      ghp=1/(1-2*j*omg2/om-omg2**2/om2)
c------------- Butterworth Tiefpass 2. Ordnung, fc=30.Hz --------------
      omg3=zpi*30.
      glp2=1/(1+1.4142*j*om/omg3-om2/omg3**2)
c----------- Butterworth Tiefpass 3.Ordnung, fc=25 Hz (Discriminator)----
      omg4=zpi*25.
      glp3=1/(1+2*j/omg4-2*om2/omg4**2-jom*om2/omg4**3)
c----------- Butterworth Tiefpass 2.Ordnung, fc=24 Hz (Antialiasing) ----
      omg5=zpi*24.
      glp2a=1/(1+1.4142*jom/omg5-om2/omg5**2)
c----------- Uebertragungsfunktion des gesamten Systems -----------------
      G=gseis*ghp*glp2*glp3*glp2a
c---- Deconvolution mit der Uebertragungsfunktion eines Wood-Anderson-Seism.
      Twa=0.8
      omw=zpi/Twa
      hw=0.78
      Vwa=2800
      G=(om2 *Vwa/(omw**2+2*jom*omw*hw-om2))/G
      ampl=cabs(G) ! [sec] Amplitudengang = Betrag der kompl. Transferfkt.
      ampl=ampl*pdrampl/10000.  ! [sec*1/10um/sec*1/1000=mm]
      return
      end
c
      subroutine ufwabo(isecpendel,seismkonst,seismdamp,voltgain,
     &                  devampl,period,iampltype,ampl,ifilter)
c---- subr. UFWABO ---------------------------------------------------
c
c     Berechnet die komplexe Uebertragungsfunktion im Frequenzbereich
c     fuer das System 'Bodenbewegung --> Devco-Auswertetisch' bzw.
c     fuer das System 'Bodenbewegung --> ADC (Mini)'
c
c     Options: iampltype = 1  --> [ampl] = mm(Wood-Anderson)/mm(Devco)
c              iampltype = 2  --> [ampl] = Mikron(Bodenbewegung)/mm(Devco)
c
c              ifilter   = 'DE' --> devampl von Develocorder
c              ifilter   = 'AD' --> devampl von ADC
c
c     Zur Subr. UFWABO gehoeren ebenfalls die Complex-Functions:
c     HP1, STS373, DEV und B6.
c--------------------------------------------------------------------------
      implicit none
      integer isecpendel,iampltype
      real seismkonst,seismdamp,voltgain,devampl
      real period,ampl
c
      real zpi,t02,t01,ts,s,hs,om,om2,omg2,tau2
      real om0,gain,twa,o0w,hw,vwa,amplresp,omg1,tau1
      integer iline,i
      character*(*) ifilter
      complex g1,g2,G,j,B6,HP1,STS373,DEV,jom
      complex gseismb
      data zpi/6.28319/j/(0.,1.)/
c
c---- Spezifikation der Hochpaesse auf der Gesamtuebertragung (PREampl. -->
c---- ADC) wie mit Programm Calap berechnet (bzw. verifiziert)
c
      iline=isecpendel
      if(iline.eq.1)then
         T02=11.9  ! 1-sec-Linie
         T01= 3.7  !     ''
      endif
      if(iline.eq.2)then
         T02= 12.7  ! 2-sec-Linie
         T01=112.   !     ''
      endif
      TS=isecpendel
      S=seismkonst
      hs=seismdamp
c
c---- Berechnung der Uebertragungsfunktion {Bodenbewegung --> ADC} ---
c
      om=zpi/period
      jom=j*om
      om2=om**2
c
c--------------------------- HP 2. Ordnung ----------------------------
c--------------------------- zus.gesetzt aus 2 HP 1. Ordnung ----------
      omg2=zpi*1./T02
      tau2=1./omg2
c
      g2=jom*tau2
      g2=(g2/(1.+g2))**2
c--------------------------- HP 1. Ordnung ----------------------------
      omg1=zpi*1./T01
      tau1=1./omg1
c
      g1=jom*tau1
      g1=g1/(1.+g1)
c----- Seismometer ----------------------------------------------------
100   om0=zpi/TS
      ! fuer Bodenbewegung ------------------------------------------
      gseismb=om2*jom*S/(om0**2+2*jom*om0*hs-om2)    ! [V/cm]
c----------------------------------------------------------------------
c
c---- Die Uebertragungsfunktionen wurden urspruenglich via ADC berechnet
c     (WK Kradi ~1983) und demzufolge ist der anti-alias-Filter schon
c     in der Uebertragungsfunktion mit inbegriffen. Um die Uebertragungs-
c     funktionen nur bis zum Diskriminator zu berechnen, muss dieser
c     Filter rueckgaengig gemacht werden:
c---- Inverse Uebertragung des DC-Unterdrueckungs-HP vor ADC
c
      G=g1*g2*gseismb/HP1(om)
c---- (damit sind wir beim Diskriminator)
c
      if(ifilter.eq.'DE')then
c
c        Uebertragung {Diskriminator --> Vergroesserungsschirm}
c
         G=G*STS373(om)*DEV(om)
                     ! Bis jetzt [mm Devco/mm Bodenbew.]
c
      endif
c
      if(ifilter.eq.'AD')then
c
c        Uebertragung {Diskriminator --> ADC }
c
         G=G*HP1(om)*B6(om)
                       ! Bis jetzt [volt/cm Bodenbewegung]
         gain=16.      ! [plot mm / volt]   (Ablesung ist so skaliert)
         G=G/10.       ! [volt/mm Bodenbewegung]
         G=G*gain      ! [mm plot/mm Bodenbewegung]
c
      endif
c
c---- Umrechnung in Mikron/Devco-Millimeter  bzw.
c---- Deconvolution mit der Uebertragungsfunktion eines Wood-Anderson-Seism.
c
      i=iampltype        !!! Ablesung = devco oder plot
      if(i.eq.2)then
         G=G/1000.     ! jetzt: mm Ablesung/um Bodenbewegung
         goto 999
      endif
      if(i.eq.1)then
         Twa=0.8
         o0w=zpi/Twa
         hw=0.78
         Vwa=2800
         G=G/(om2 *Vwa/(o0w**2+2*jom*o0w*hw-om2))
C                                             ! [1]     nur Wegvergr. !!
                     ! jetzt: mm Ablesung/mm Wood-Anderson !
         goto 999
      endif
c
999   amplresp=cabs(G)   ! Amplitudengang := Betrag der kompl. Transferfkt.
c
      ampl=devampl/(voltgain*amplresp)     ! erst hier reziproke Einheiten !!
c
c fuer IAMPLTYPE =          1          |           2
c dann [ampl]    = [mm Wood-Anderson]  |  [um Bodenbewegung]
c
      RETURN
      end
c
c
c------------------- B6: Bessel-Filter (TP) ------------------------------
      complex function B6(om)
      implicit none
      complex j,gbn1,gbn2,gbn3
      real om,zpi,fg1,omg1
c
      data j/(0.,1.)/,zpi/6.28319/
      fg1=12.     ! Grenzfrequenz 12 Hz
      omg1=om/(zpi*fg1)
c
      gbn1= 1. + 1.2217*j*omg1 - 0.3887*omg1*omg1
      gbn2= 1. + 0.9686*j*omg1 - 0.3505*omg1*omg1
      gbn3= 1. + 0.5131*j*omg1 - 0.2756*omg1*omg1
c
      b6=1./(gbn1*gbn2*gbn3)
c
      end
c
c------------------- HP 0.1 Hz 1. Ordnung (vor ADC) [DC-Unterdrueckung]-
      complex function HP1(om)
      implicit none
      real om
      complex j
      real zpi,fg,tau
      data j/(0.,1.)/zpi/6.28319/
      Fg=0.1 ! [Hz]     Grenzfrequenz
      tau=1./(zpi*Fg)
c
      HP1=j*om*tau
      HP1=HP1/(1.+HP1)
c
      return
      end
c
c------------------- TP STS-373 (Bessel 4. Ordnung, 5 Hz ) ----------
      complex function STS373(om)
      implicit none
      real om
c      implicit InTeGeR*4 (i-n)
      complex    j,gbn1,gbn2, o
      real zpi,fg1,omg1
      data j/(0.,1.)/zpi/6.28319/
      fg1=5.
      omg1=zpi*fg1
      o=om/omg1
c
      gbn1=1.+1.3397*j*o-0.4889*o**2
      gbn2=1.+0.7743*j*o-0.3890*o**2
c
      STS373=1/(gbn1*gbn2)
c
      return
      end
c
c------------------- Develocorder & 12-fache Schirmvergroesserung -----
      complex function DEV(om)
      implicit none
      real om
      complex j
      real zpi,fg,tau
      data j/(0.,1.)/zpi/6.28319/
                                 ! 1.6 cm/V MIT 12-f. Vergr. am Tisch !!
      Fg=15. ! [Hz]              ! Galvanometer      krit. Daempfung
      tau=1./(zpi*Fg)            ! (TP)              2. Ordnung
c
      DEV=1.6*( 1./(1.+j*om*tau) )**2           ! [cm/V]
c
c!!!!!!!! bei neuen functions diese in Haupt-Subroutine deklarieren  !!!
      return
      end
c
cek   end of vel_mag.f
c
cek   begin of vel_math.f
c
      subroutine STOREG(k,l)
c
c---- save column k of g in gcopy
c
      implicit none
      integer k,l
      include 'vel_com.f'
      integer ii,jj,i,j
c
      ii=0
      jj=0
      do 1 i=1,nvar
      do 1 j=1,i
      jj=jj+1
      if(i.ne.k.and.j.ne.k) goto 1
      ii=ii+1
      gcopy(l,ii)=g(jj)
   1  continue
      return
      end ! of subr. storeg
c
      subroutine SINGULARVALUES(i)
c
c     compute the singular values (in this case: eigenvalues) of the
c     symmetric matrix G.
c
      implicit none
      integer i
      include 'vel_com.f'
      integer k,ii,jj,i1,j
c
      real A(4,4), Xsol(4)
c
      do k=1,nvar
         ii=0
         jj=0
         do i1=1,nvar
            do j=1,i1
               jj=jj+1
               if(.NOT.(i1.ne.k.and.j.ne.k))then
                  ii=ii+1
                  A(k,ii)=g(jj)
               endif
            enddo
         enddo
      enddo
c
      call SVDSOLUK(A,RHt,nvar,-1.,Xsol,SV,ale(i),COVs,Rs)
c
      if(ifixsolution.eq.1) call ALESUBR(SV,3,ale(i))
      if(ifixsolution.eq.9) call ALESUBR(SV,1,ale(i))
c
      if(.not.single_turbo)then
         write(16,*)'Singular values; iteration #',nitt
         write(16,'(1x,4(2x,f10.6))') (sv(jj),jj=1,nvar)
         write(16,*)'ALE = ',ale(i)
      endif
c
ccc      write(16,*)'SVD-solution vector:'
ccc      write(16,'(1x,4(2x,f16.6))') (Xsol(jj),jj=1,nvar)
c
      return
      end ! of subr. singularvalues
c
      subroutine DAMPG
c
c     apply damping to the diagonal elements of the symmetric matrix G
c
      implicit none
      include 'vel_com.f'
      integer j,k
c
c     apply damping to diagonal elements of G-matrix :
c
      j=0
      do k=1,nvar          !   k = 1 2 3  4  5 ...
         j=j+k             !   j = 1 3 6 10 15 ...   <== diagonal elements of G
         g(j)=g(j)+othet
      enddo
c
      return
      end ! of subr. dampg
c
      subroutine ACCUNORMEQS(rowofa,nvar,res,w,g,rhs)
c
c     Urs Kradolfer, 22. 4. 1987
c
c     accumulate normal equations, knowing only one row of A at a time;
c     subr. ACCUNORMEQS must be called nvar times (once for
c     each observation).
c     w is the (normalized ! ) weight of the observation.
c
c        A * x = res          INPUT: one row of A and one element of res
c
c     At*A * x = At*res
c
c       G  * x =  RHS         OUTPUT: G (nvar x nvar) and RHS (nvar)
c
      implicit none
      real g(*),rowofa(*),rhs(*),res,w
      integer nvar
      call OUTER(g,rowofa,nvar,w)
      call RSIDE(rhs,rowofa,nvar,res,w)
c
      return
      end ! of subr. accunormeqs
c
      subroutine OUTER(g,s,n,w)
c
c   this routine does an outer product of a vector
c   s with itself.  This is used in least squares
c   to accumulate (At)A knowing only one
c   row of A at a time.  (At)A is a symmetric
c   matrix which is n by n.  This matrix is stored
c   in the vector g in symmetric storage mode.  Thus,
c   given the ith row and jth column of the matrix, then
c            nsym=((i*(i-1))/2) + j
c   gives the index for the element nsym of g
c   for i.ge.j.  For j less than i, the
c   element ij is identical to the element ji.
c
c   input:
c    s(i) one row of the matrix (A) .
c    n - the length of s(i)
c    w - a weight for the least squares.
c
      implicit none
      integer n
      real s(n),g((n*(n-1)/2)+n),w
      integer i,j,nsym
      real a,b
c
      do 1 i=1,n
      if(s(i).eq.0) goto 1
      a=s(i)*w
c
      do 2 j=1,i
      if(s(j).eq.0) goto 2
      b=s(j)*a
c
      nsym=((i*(i-1))/2) + j
      g(nsym)=g(nsym) + b
  2   continue
c
  1   continue
      return
      end ! of subr. outer
c
      subroutine RSIDE(rht,s,n,res,w)
c
c   this routine accumulates in rht the quantity (At)b
c   for the normal equations.
c
c   input:
c    s(i)-one column of the matrix (At)
c    n - the length of s(i)
c    res - one value of the vector b.
c    w  -  a weight for the weighted least sqs.
c
      implicit none
      integer n
      real rht(n),s(n),w,res
      integer i
c
      do 1 i=1,n
      if(s(i).eq.0) goto 1
      rht(i)=res*w*s(i) + rht(i)
1     continue
      return
      end ! of subr. rside
c
      subroutine ALEsubr(SV,m,ale)
c
c     Urs Kradolfer, June 1st, 1987
c
c     Input : Array SV, containing the M singular values
c             m , dimension of SV
c
c     Output: ALE-value
c
c             ALE = - 1./(m-izero) * sum(i=1...m) log10( SV(i)/SV(1) )
c                   + izero*10.
c
c             where izero = # zero-singular-values
c
      implicit none
c
      integer m
      real SV(m)
      real ale,alesum
      integer i,izero
c
      if(m.eq.0)stop'ALEsubr>>> m was zero !!!'
c
      ale=0.0
      alesum=0.0
      izero=0
      do i=1,m
         ale=sv(i)/sv(1)
         if(ale.le.0.0)then
            izero=izero+1
         else
            alesum=alesum+log10(ale)
         endif
      enddo
      ale=-(alesum/float(m-izero)) + 10*izero
c
      return
      end ! of subr. alesubr
c
      subroutine LSVDF  (A,IA,M,N,B,IB,NB,S,WK,IER)
c     This is an IMSL-subroutine
      implicit none
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,IB,NB,IER
      real               A(IA,N),B(IB,1),S(N),WK(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JP1,K,L,MM,NN,NNP1,NS,NSP1
      real               ZERO,ONE,T
      integer            nm
      real               f
      DATA               ZERO/0.0e0/,ONE/1.0e0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
C                                  BEGIN SPECIAL FOR ZERO ROWS AND
C                                    COLS. PACK THE NONZERO COLS TO THE
C                                    LEFT
      NN=N
      IER=34
      IF (NN.LE.0.OR.M.LE.0) GO TO 9000
      IER=0
      J=NN
    5 CONTINUE
      DO 10 I=1,M
         IF (A(I,J).NE.ZERO) GO TO 25
   10 CONTINUE
C                                  COL J IS ZERO. EXCHANGE IT WITH COL
C                                    N
      IF (J.EQ.NN) GO TO 20
      DO 15 I=1,M
   15 A(I,J)=A(I,NN)
   20 CONTINUE
      A(1,NN)=J
      NN=NN-1
   25 CONTINUE
      J=J-1
      IF (J.GE.1) GO TO 5
C                                  IF N=0 THEN A IS ENTIRELY ZERO AND
C                                    SVD COMPUTATION CAN BE SKIPPED
      NS=0
      IF (NN.EQ.0) GO TO 120
C                                  PACK NONZERO ROWS TO THE TOP QUIT
C                                    PACKING IF FIND N NONZERO ROWS
      I=1
      MM=M
   30 IF (I.GT.N.OR.I.GE.MM) GO TO 75
      IF (A(I,I).NE.ZERO) GO TO 40
      DO 35 J=1,NN
         IF (A(I,J).NE.ZERO) GO TO 40
   35 CONTINUE
      GO TO 45
   40 I=I+1
      GO TO 30
C                                  ROW I IS ZERO EXCHANGE ROWS I AND M
   45 IF (NB.LE.0) GO TO 55
      DO 50 J=1,NB
         T=B(I,J)
         B(I,J)=B(MM,J)
         B(MM,J)=T
   50 CONTINUE
   55 DO 60 J=1,NN
   60 A(I,J)=A(MM,J)
      IF (MM.GT.NN) GO TO 70
      DO 65 J=1,NN
   65 A(MM,J)=ZERO
   70 CONTINUE
C                                  EXCHANGE IS FINISHED
      MM=MM-1
      GO TO 30
C
   75 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
C                                  BEGIN SVD ALGORITHM..
C                                  (1) REDUCE THE MATRIX TO UPPER
C                                    BIDIAGONAL FORM WITH HOUSEHOLDER
C                                    TRANSFORMATIONS.
C                                    H(N)...H(1)AQ(1)...Q(N-2) =
C                                    (D**T,0)**T WHERE D IS UPPER
C                                    BIDIAGONAL.
C                                  (2) APPLY H(N)...H(1) TO B. HERE
C                                    H(N)...H(1)*B REPLACES B IN
C                                    STORAGE.
C                                  (3) THE MATRIX PRODUCT W=
C                                    Q(1)...Q(N-2) OVERWRITES THE FIRST
C                                    N ROWS OF A IN STORAGE.
C                                  (4) AN SVD FOR D IS COMPUTED. HERE K
C                                    ROTATIONS RI AND PI ARE COMPUTED
C                                    SO THAT RK...R1*D*P1**(T)...PK**(T)
C                                    = DIAG(S1,...,SM) TO WORKING
C                                    ACCURACY. THE SI ARE NONNEGATIVE
C                                    AND NONINCREASING. HERE RK...R1*B
C                                    OVERWRITES B IN STORAGE WHILE
C                                    A*P1**(T)...PK**(T) OVERWRITES A
C                                    IN STORAGE.
C                                  (5) IT FOLLOWS THAT,WITH THE PROPER
C                                    DEFINITIONS, U**(T)*B OVERWRITES
C                                    B, WHILE V OVERWRITES THE FIRST N
C                                    ROW AND COLUMNS OF A.
      L=MIN0(MM,NN)
C                                  THE FOLLOWING LOOP REDUCES A TO
C                                    UPPER BIDIAGONAL AND ALSO APPLIES
C                                    THE PREMULTIPLYING TRANSFORMATIONS
C                                    TO B.
      DO 85 J=1,L
         IF (J.GE.MM) GO TO 80
         JP1=MIN0(J+1,NN)
         CALL VHS12 (1,J,J+1,MM,A(1,J),1,T,A(1,JP1),1,IA,NN-J)
         CALL VHS12 (2,J,J+1,MM,A(1,J),1,T,B,1,IB,NB)
   80    IF (J.GE.NN-1) GO TO 85
         CALL VHS12 (1,J+1,J+2,NN,A(J,1),IA,WK(J,2),A(J+1,1),IA,1,MM-J)
   85 CONTINUE
C                                  COPY THE BIDIAGONAL MATRIX INTO THE
C                                    ARRAY S FOR LSVDB
      IF (L.EQ.1) GO TO 95
      DO 90 J=2,L
         S(J)=A(J,J)
         WK(J,1)=A(J-1,J)
   90 CONTINUE
   95 S(1)=A(1,1)
C
      NS=NN
      IF (MM.GE.NN) GO TO 100
      NS=MM+1
      S(NS)=ZERO
      WK(NS,1)=A(MM,MM+1)
  100 CONTINUE
C                                  CONSTRUCT THE EXPLICIT N BY N
C                                    PRODUCT MATRIX, W=Q1*Q2*...*QL*I
C                                    IN THE ARRAY A
      DO 115 K=1,NN
         I=NN+1-K
         IF (I.GT.MIN0(MM,NN-2)) GO TO 105
         CALL VHS12 (2,I+1,I+2,NN,A(I,1),IA,WK(I,2),A(1,I+1),1,IA,NN-I)
  105    DO 110 J=1,NN
  110    A(I,J)=ZERO
         A(I,I)=ONE
  115 CONTINUE
C                                  COMPUTE THE SVD OF THE BIDIAGONAL
C                                    MATRIX
C
C      LEVEL=1
C      CALL UERSET(LEVEL,LEVOLD)
      CALL LSVDB (S(1),WK(1,1),NS,A,IA,NN,B,IB,NB,IER)
C                                  TEST FOR IER=33
C
      IF (IER.GT.128) GO TO 9000
C      CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.NE.33) GO TO 120
cuk      WRITE(6,4500) IER
 4500 FORMAT(2X,'SUBR. LSVDF: IER=',I4,'AFTER CALL OF SUBR.LSVDB')
      T=0.0e0
      NM=MIN0(M,N)
      IF (S(1).NE.ZERO) T=S(NM)/S(1)
      F=100.0e0+T
      IF (F.EQ.100.0e0) GO TO 120
      IER=0
  120 CONTINUE
      IF (NS.GE.NN) GO TO 130
      NSP1=NS+1
      DO 125 J=NSP1,NN
  125 S(J)=ZERO
  130 CONTINUE
      IF (NN.EQ.N) GO TO 155
      NNP1=NN+1
C                                  MOVE RECORD OF PERMUTATIONS AND
C                                    STORE ZEROS
      DO 140 J=NNP1,N
         S(J)=A(1,J)
         IF (NN.LT.1) GO TO 140
         DO 135 I=1,NN
  135    A(I,J)=ZERO
  140 CONTINUE
C                                  PERMUTE ROWS AND SET ZERO SINGULAR
C                                    VALUES
      DO 150 K=NNP1,N
         I=S(K)
         S(K)=ZERO
         DO 145 J=1,N
            A(K,J)=A(I,J)
  145    A(I,J)=ZERO
         A(I,K)=ONE
  150 CONTINUE
C                                  END SPECIAL FOR ZERO ROWS AND
C                                    COLUMNS
  155 IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
C      CALL UERTST (IER,'LSVDF ')
cuk
      if(ier.ne.33)then
      WRITE(6,4501) IER
 4501 FORMAT(2X,'SUBR. LSVDF: AT END, IER=',I4,'CHECK SUBR. HEAD FOR',
     11X,'MEANING')
      endif
cuk
 9005 RETURN
      END
c
      subroutine LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
c
C   IMSL ROUTINE NAME   - LSVDB
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - HP1000/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SINGULAR VALUE DECOMPOSITION OF A BIDIAGONAL
C                           MATRIX.
C
C   USAGE               - CALL LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C
C   ARGUMENTS    D      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                         ON INPUT, D CONTAINS THE DIAGONAL ELEMENTS
C                           OF THE BIDIAGONAL MATRIX B. D(I)=B(I,I),
C                           I=1,...,N.
C                         ON OUTPUT, D CONTAINS THE N (NONNEGATIVE)
C                           SINGULAR VALUES OF B IN NONINCREASING
C                           ORDER.
C                E      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                         ON INPUT, E CONTAINS THE SUPERDIAGONAL
C                           ELEMENTS OF B. E(1) IS ARBITRARY,
C                           E(I)=B(I-1,I), I=2,...,N.
C                         ON OUTPUT, THE CONTENTS OF E ARE MODIFIED
C                           BY THE SUBROUTINE.
C                N      - ORDER OF THE MATRIX B. (INPUT)
C                V      - NRV BY N MATRIX. (INPUT/OUTPUT)
C                           IF NRV.LE.0, V IS NOT USED. OTHERWISE,
C                           V IS REPLACED BY THE NRV BY N PRODUCT
C                           MATRIX V*VB. SEE REMARKS.
C                IV     - ROW DIMENSION OF MATRIX V EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NRV    - NUMBER OF ROWS OF V. (INPUT)
C                C      - N BY NCC MATRIX. (INPUT/OUTPUT)
C                           IF NCC.LE.0 C IS NOT USED. OTHERWISE, C
C                           IS REPLACED BY THE N BY NCC PRODUCT
C                           MATRIX UB**(T) * C. SEE REMARKS.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NCC    - NUMBER OF COLUMNS IN C. (INPUT)
C                IER    - ERROR PARAMETER. (INPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT MATRIX B IS NOT FULL
C                             RANK OR VERY ILL-CONDITIONED. SMALL
C                             SINGULAR VALUES MAY NOT BE VERY ACCURATE.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CONVERGENCE WAS
C                             NOT ATTAINED AFTER 10*N QR SWEEPS.
C                             (CONVERGENCE USUALLY OCCURS IN ABOUT
C                             2*N SWEEPS).
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - LSVG1,LSVG2,VHS12,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      LSVDB COMPUTES THE SINGULAR VALUE DECOMPOSITION OF
C                AN N BY N BIDIAGONAL MATRIX
C                     B = UB * S * VB**(T)    WHERE
C                UB AND VB ARE N BY N ORTHOGONAL MATRICES AND
C                S IS DIAGONAL.
C                IF ARGUMENTS V AND C ARE N BY N IDENTITY MATRICES,
C                ON EXIT THEY ARE REPLACED BY VB AND UB**T,
C                RESPECTIVELY.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    N,IV,NRV,IC,NCC,IER
cuk      real    D(N),E(N),V(IV,1),C(IC,1)
      real       D(N),E(N),V(IV,100),C(IC,100)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I,II,J,K,KK,L,LL,LP1,NQRS,N10
      LOGICAL    WNTV,HAVERS,FAIL
      real       DNORM,ZERO,ONE,TWO,CS,F,G,H,SN,T,X,Y,Z
      DATA       ZERO/0.0e0/,ONE/1.0e0/,TWO/2.0e0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (N.LE.0) GO TO 9005
      N10=10*N
      WNTV=NRV.GT.0
      HAVERS=NCC.GT.0
      FAIL=.FALSE.
      NQRS=0
      E(1)=ZERO
      DNORM=ZERO
      DO 5 J=1,N
    5 DNORM=MAX1(ABS(D(J))+ABS(E(J)),DNORM)
      DO 100 KK=1,N
         K=N+1-KK
C                                  TEST FOR SPLITTING OR RANK
C                                    DEFICIENCIES FIRST MAKE TEST FOR
C                                    LAST DIAGONAL TERM, D(K), BEING
C                                    SMALL.
   10    IF (K.EQ.1) GO TO 25
         T=DNORM+D(K)
         IF (T.NE.DNORM) GO TO 25
C
C                                  SINCE D(K) IS SMALL WE WILL MAKE A
C                                    SPECIAL PASS TO TRANSFORM E(K) TO
C                                    ZERO.
         CS=ZERO
         SN=-ONE
         DO 20 II=2,K
            I=K+1-II
            F=-SN*E(I+1)
            E(I+1)=CS*E(I+1)
            T=D(I)
            CALL LSVG1 (T,F,CS,SN,D(I))
C                                  TRANSFORMATION CONSTRUCTED TO ZERO
C                                    POSITION (I,K).
            IF (.NOT.WNTV) GO TO 20
            DO 15 J=1,NRV
   15       CALL LSVG2 (CS,SN,V(J,I),V(J,K))
C
C                                  ACCUMULATE RT. TRANSFORMATIONS IN V.
   20    CONTINUE
C                                  THE MATRIX IS NOW BIDIAGONAL, AND OF
C                                    LOWER ORDER SINCE E(K) .EQ. ZERO
   25    DO 30 LL=1,K
            L=K+1-LL
            T=DNORM+E(L)
            IF (T.EQ.DNORM) GO TO 50
            T=DNORM+D(L-1)
            IF (T.EQ.DNORM) GO TO 35
   30    CONTINUE
C                                  THIS LOOP CANT COMPLETE SINCE E(1) =
C                                    ZERO.
         GO TO 50
C                                  CANCELLATION OF E(L), L.GT.1.
   35    CS=ZERO
         SN=-ONE
         DO 45 I=L,K
            F=-SN*E(I)
            E(I)=CS*E(I)
            T=DNORM+F
            IF (T.EQ.DNORM) GO TO 50
            T=D(I)
            CALL LSVG1 (T,F,CS,SN,D(I))
            IF (.NOT.HAVERS) GO TO 45
            DO 40 J=1,NCC
   40       CALL LSVG2 (CS,SN,C(I,J),C(L-1,J))
   45    CONTINUE
C                                  TEST FOR CONVERGENCE
   50    Z=D(K)
         IF (L.EQ.K) GO TO 85
C                                  SHIFT FROM BOTTOM 2 BY 2 MINOR OF
C                                    B**(T)*B.
         X=D(L)
         Y=D(K-1)
         G=E(K-1)
         H=E(K)
         F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
         G=SQRT(ONE+F**2)
         IF (F.LT.ZERO) GO TO 55
         T=F+G
         GO TO 60
   55    T=F-G
   60    F=((X-Z)*(X+Z)+H*(Y/T-H))/X
C                                  NEXT QR SWEEP
         CS=ONE
         SN=ONE
         LP1=L+1
         DO 80 I=LP1,K
            G=E(I)
            Y=D(I)
            H=SN*G
            G=CS*G
            CALL LSVG1 (F,H,CS,SN,E(I-1))
            F=X*CS+G*SN
            G=-X*SN+G*CS
            H=Y*SN
            Y=Y*CS
            IF (.NOT.WNTV) GO TO 70
C                                  ACCUMULATE ROTATIONS (FROM THE
C                                    RIGHT) IN V
            DO 65 J=1,NRV
   65       CALL LSVG2 (CS,SN,V(J,I-1),V(J,I))
   70       CALL LSVG1 (F,H,CS,SN,D(I-1))
            F=CS*G+SN*Y
            X=-SN*G+CS*Y
            IF (.NOT.HAVERS) GO TO 80
            DO 75 J=1,NCC
   75       CALL LSVG2 (CS,SN,C(I-1,J),C(I,J))
C
C                                  APPLY ROTATIONS FROM THE LEFT TO
C                                    RIGHT HAND SIDES IN C
   80    CONTINUE
         E(L)=ZERO
         E(K)=F
         D(K)=X
         NQRS=NQRS+1
         IF (NQRS.LE.N10) GO TO 10
C                                  RETURN TO TEST FOR SPLITTING.
         FAIL=.TRUE.
C                                  CUTOFF FOR CONVERGENCE FAILURE. NQRS
C                                    WILL BE 2*N USUALLY.
   85    IF (Z.GE.ZERO) GO TO 95
         D(K)=-Z
         IF (.NOT.WNTV) GO TO 95
         DO 90 J=1,NRV
   90    V(J,K)=-V(J,K)
   95    CONTINUE
C                                  CONVERGENCE. D(K) IS MADE
C                                    NONNEGATIVE
  100 CONTINUE
      IF (N.EQ.1) GO TO 140
      DO 105 I=2,N
         IF (D(I).GT.D(I-1)) GO TO 110
  105 CONTINUE
      GO TO 140
C                                  EVERY SINGULAR VALUE IS IN ORDER
  110 DO 135 I=2,N
         T=D(I-1)
         K=I-1
         DO 115 J=I,N
            IF (T.GE.D(J)) GO TO 115
            T=D(J)
            K=J
  115    CONTINUE
         IF (K.EQ.I-1) GO TO 135
         D(K)=D(I-1)
         D(I-1)=T
         IF (.NOT.HAVERS) GO TO 125
         DO 120 J=1,NCC
            T=C(I-1,J)
            C(I-1,J)=C(K,J)
  120    C(K,J)=T
  125    IF (.NOT.WNTV) GO TO 135
         DO 130 J=1,NRV
            T=V(J,I-1)
            V(J,I-1)=V(J,K)
  130    V(J,K)=T
  135 CONTINUE
C                                  END OF ORDERING ALGORITHM.
  140 IER=129
      IF (FAIL) GO TO 9000
C                                  CHECK FOR POSSIBLE RANK DEFICIENCY
      IER=33
      T=0.0e0
      IF (D(1).NE.ZERO) T=D(N)/D(1)
      F=100.0e0+T
      IF (F.EQ.100.0e0) GO TO 9000
      IER=0
      GO TO 9005
 9000 CONTINUE
C      CALL UERTST (IER,'LSVDB ')
cuk
      if(ier.ne.33)then
      WRITE(6,4500) IER
 4500 FORMAT(2X,'SUBR. LSVDB AT END: IER=',I4,' CHECK HEADER FOR',
     1' MEANING')
      endif
cuk
 9005 RETURN
      END
c
      subroutine LSVG1  (A,B,DCOS,DSIN,SIG)
c
C   IMSL ROUTINE NAME   - LSVG1
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - HP1000/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      real                         A,B,DCOS,DSIN,SIG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      real                         AA,BB
C                                  FIRST EXECUTABLE STATEMENT
      double precision doubleb, doubleaa, doubler
      IF (ABS(A).LE.ABS(B)) GO TO 5
      AA=ABS(A+A)
cuk   next statement caused an underflow! u.kradolfer, 16.7.81
cuk      SIG=AA*SQRT(0.25e0+(B/AA)**2)
cuk   next statement produced a compiler warning with range checking option
cuk      SIG=AA*sngl( dSQRT(0.25d0+(dble(B)/dble(AA))**2) ) ! u.k. 7.2.92
      doubleb=DBLE(b)
      doubleaa=DBLE(aa)
      doubler=dSQRT(0.25d0+(doubleb/doubleaa)**2)
      sig=AA*sngl(doubler)
      DCOS=A/SIG
      DSIN=B/SIG
      RETURN
    5 IF (B.EQ.0.0e0) GO TO 10
      BB=ABS(B+B)
      SIG=BB*SQRT(0.25e0+(A/BB)**2)
      DCOS=A/SIG
      DSIN=B/SIG
      RETURN
   10 SIG=0.0e0
      DCOS=0.0e0
      DSIN=1.0e0
      RETURN
      END
c
      subroutine LSVG2  (DCOS,DSIN,X,Y)
c
C   IMSL ROUTINE NAME   - LSVG2
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - HP1000/DOUBLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      real                         DCOS,DSIN,X,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      real                         XR
C                                  FIRST EXECUTABLE STATEMENT
      XR=DCOS*X+DSIN*Y
      Y=-DSIN*X+DCOS*Y
      X=XR
      RETURN
      END
c
      subroutine VHS12  (MODE,LP,L1,M,U,INCU,UP,C,INCC,ICV,NCV)
C
      implicit none
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MODE,LP,L1,M,INCU,INCC,ICV,NCV
cuk      real               U(1),UP,C(1)
cuk      real               U(m+1),UP,C(ncv*m)  ! did NOT work...
      real               U(100),UP,C(100)   ! in this case only!!! (single-ev-l)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IJ,ILP,IL1,IM,INCR,I2,I3,I4,J
      real               SM,B
      real               ONE,CL,CLINV,SM1
C                                  FIRST EXECUTABLE STATEMENT
      ONE=1.e0
C
      IF (0.GE.LP.OR.LP.GE.L1.OR.L1.GT.M) GO TO 9005
      ILP=(LP-1)*INCU+1
      IL1=(L1-1)*INCU+1
      IM=(M-1)*INCU+1
      CL=ABS(U(ILP))
      IF (MODE.EQ.2) GO TO 15
C                                  CONSTRUCT THE TRANSFORMATION.
      DO 5 IJ=IL1,IM,INCU
    5 CL=MAX1(ABS(U(IJ)),CL)
      IF (CL.LE.0.0e0) GO TO 9005
      CLINV=ONE/CL
      SM=(U(ILP)*CLINV)**2
      DO 10 IJ=IL1,IM,INCU
   10 SM=SM+(U(IJ)*CLINV)**2
C                                  CONVERT DBLE. PREC. SM TO SNGL.
C                                    PREC. SM1
      SM1=SM
      CL=CL*SQRT(SM1)
      IF (U(ILP).GT.0.0e0) CL=-CL
      UP=U(ILP)-CL
      U(ILP)=CL
      GO TO 20
C                                  APPLY THE TRANSFORMATION
C                                    I+U*(U**T)/B TO C.
   15 IF (CL.LE.0.0e0) GO TO 9005
   20 IF (NCV.LE.0) GO TO 9005
      B=UP*U(ILP)
C                                  B MUST BE NONPOSITIVE HERE. IF B =
C                                    0., RETURN.
      IF (B.GE.0.0e0) GO TO 9005
      B=ONE/B
      I2=1-ICV+INCC*(LP-1)
      INCR=INCC*(L1-LP)
      DO 35 J=1,NCV
         I2=I2+ICV
         I3=I2+INCR
         I4=I3
         SM=C(I2)*UP
         DO 25 IJ=IL1,IM,INCU
            SM=SM+C(I3)*U(IJ)
            I3=I3+INCC
   25    CONTINUE
         IF (SM.EQ.0.0e0) GO TO 35
         SM=SM*B
         C(I2)=C(I2)+SM*UP
         DO 30 IJ=IL1,IM,INCU
            C(I4)=C(I4)+SM*U(IJ)
            I4=I4+INCC
   30    CONTINUE
   35 CONTINUE
 9005 RETURN
      END
c
      subroutine LUDECP (a,ul,n,d1,d2,ier)
c
c     function:        - cholesky decomposition of a matrix A -
c     ---------
c     A   nxn pos. def. symm. matrix, stored in symm. storage mode
c     UL  nxn lower triangular matrix,  with UL * ULt = A
c     !!! diagonal elements of UL are stored in reziprocal form !!!
c     det(A) = D1*(2.**D2)      d1, d2 is output of LUDECP
c     IER=0  --> A is pos. def.;    everything o.k.
c
      implicit none
      integer n
      integer ier
      real   d1,d2
      real   a(n*(n+1)/2),ul(n*(n+1)/2)
      real   zero,one,four,sixtn,sixth,rn,x
      integer ip,i,iq,ir,j,k,ip1
      data   zero,one,four,sixtn,sixth/0.0,1.,4.,16.,.0625/
c
      d1=one
      d2=zero
      rn=one/(n*sixtn)
      ip=1
      ier=0
      do 45 i=1,n
         iq=ip
         ir=1
         do 40 j=1,i
            x=a(ip)
            if (j .eq. 1) go to 10
            do 5  k=iq,ip1
               x=x-ul(k)*ul(ir)
               ir=ir+1
    5       continue
   10       if (i.ne.j) go to 30
            d1=d1*x
            if (a(ip)+x*rn .le. a(ip)) go to 50
   15       if (abs(d1) .le. one) go to 20
            d1=d1 * sixth
            d2=d2 + four
            go to 15
   20       if (abs(d1) .ge. sixth) go to 25
            d1=d1 * sixtn
            d2=d2 - four
            go to 20
   25       ul(ip)=one/sqrt(x)
            go to 35
   30       ul(ip)=x * ul(ir)
   35       ip1=ip
            ip=ip+1
            ir=ir+1
   40    continue
   45 continue
      go to 9005
   50 ier=129
 9000 continue
 9005 return
      end ! of subr. ludecp
c
      subroutine luelmp (a,b,n,x)
c
c     function:     - elimination part of the solution of Ax=b
c     ---------
c       A  (nxn) lower triangular matrix (outout of LUDECP)
c       !!! diagonal elements of A are stored in reziprocal form !!!
c       b  n-vector
c       x  n-vector (solution of equation Ax=b)
c
      implicit none
      integer n
      real a(n*(n+1)/2),b(n),x(n)
      real zero,t
      integer ip,iw,i,im1,k,n1,ii,is,iq,kk
      data zero/0./
c                              solution of ly = b
      ip=1
      iw=0
      do 15 i=1,n
         t=b(i)
         im1=i-1
         if (iw .eq. 0) go to 9
         ip=ip+iw-1
         do 5 k=iw,im1
            t=t-a(ip)*x(k)
            ip=ip+1
    5    continue
         go to 10
    9    if (t .ne. zero) iw=i
         ip=ip+im1
   10    x(i)=t*a(ip)
         ip=ip+1
   15 continue
c                                  solution of ux = y
      n1=n+1
      do 30 i=1,n
         ii=n1-i
         ip=ip-1
         is=ip
         iq=ii+1
         t=x(ii)
         if (n.lt.iq) go to 25
         kk=n
         do 20 k=iq,n
            t=t-a(is)*x(kk)
            kk=kk-1
            is=is-kk
   20    continue
   25    x(ii)=t*a(is)
   30 continue
      return
      end ! of subr. luelmp
c
      subroutine MATRMULT(A,m,p,B,p1,n,C,m1,n1)
c     Urs Kradolfer, 28.3.1987
      implicit none
c     input:  (mxp)-matrix A
c             (p1xn)-matrix B
c     output: (m1xn1)-matrix C = A*B
      integer m,p, p1,n, m1,n1 ,
     &        i,j,k
      real A(m,p), B(p1,n), C(m1,n1) ,
     &     s
      if(m.ne.m1)stop'Matrices cannot be multiplied !'
      if(p.ne.p1)stop'Matrices cannot be multiplied !'
      if(n.ne.n1)stop'Matrices cannot be multiplied !'
      do i=1,m
         do j=1,n
            s=0.
            do k=1,p
               s=s+A(i,k)*B(k,j)
            enddo
            C(i,j)=s
         enddo
      enddo
      return
      end ! of subr. matrmult
c
      subroutine MATRTRAN(A,n,m,AT)
c     Urs Kradolfer, 28.3.1987
      implicit none
c     input:  (nxm)-matrix A
c     output: (mxn)-matrix AT ( = A transpose )
      integer n,m ,
     &        i,j
      real A(n,m),AT(m,n)
      do i=1,n
         do j=1,m
            AT(j,i)=A(i,j)
         enddo
      enddo
      return
      end ! of subr. matrtran
c
      subroutine MAXII(n,nx,imax,jndex)
c determine maximum-value of an integer-array
      implicit none
      integer jndex,i,n,imax, nx(n)
      jndex=1
      do 1 i=1,n
   1  if(nx(jndex).le.nx(i)) jndex=i
      imax=nx(jndex)
      return
      end ! of subr. maxni
c
      subroutine MAXRI(n,x,xmax,jndex)
c determine maximum-value of a real-array
      implicit none
      integer i,n,jndex
      real x(n), xmax
      jndex=1
      do 1 i=1,n
   1  if(x(jndex).le.x(i)) jndex=i
      xmax=x(jndex)
      return
      end ! of subr. maxri
c
      subroutine SORTI(iX,NO)   ! Urs Kradolfer, Okt. 1986 (aus HYPOSUBR.FOR )
c
c  SORTA sortiert einen Integer-Array iX so, dass iX(1) < iX(2) < ... < iX(N)
c
c  Aufruf: call SORTI(iX,N)
c
c  mit INTEGER iX(N)
c      INTEGER N       (Anzahl Elemente in X)
c
c  Array iX wird in SORTI veraendert !!
c
      implicit none
      integer i,no,mo,ko,jo
      integer iX(no),itemp
c
      if(no.eq.0) return
      MO=NO
 2    IF (MO-15) 21,21,23
 21   IF (MO-1) 29,29,22
 22   MO=2*(MO/4)+1
      GO TO 24
 23   MO=2*(MO/8)+1
 24   KO=NO-MO
      JO=1
 25   I=JO
 26   IF (iX(I)-iX(I+MO)) 28,28,27
 27   iTEMP=iX(I)
      iX(I)=iX(I+MO)
      iX(I+MO)=iTEMP
      I=I-MO
      IF (I-1) 28,26,26
 28   JO=JO+1
      IF (JO-KO) 25,25,2
 29   continue
      return
      END ! of subr. sorti
c
      subroutine SORTR(X,NO)   ! Urs Kradolfer, Okt. 1986 (aus HYPOSUBR.FOR )
c
c  SORTR sortiert einen Real-Array X so, dass X(1) < X(2) < ... < X(N)
c
c  Aufruf: call SORTR(X,N)
c
c  mit REAL    X(N)
c      INTEGER N       (Anzahl Elemente in X)
c
c  Array X wird in SORTR veraendert !!
c
      implicit none
      integer i,no,mo,ko,jo
      real X(no),temp
c
      MO=NO
 2    IF (MO-15) 21,21,23
 21   IF (MO-1) 29,29,22
 22   MO=2*(MO/4)+1
      GO TO 24
 23   MO=2*(MO/8)+1
 24   KO=NO-MO
      JO=1
 25   I=JO
 26   IF (X(I)-X(I+MO)) 28,28,27
 27   TEMP=X(I)
      X(I)=X(I+MO)
      X(I+MO)=TEMP
      I=I-MO
      IF (I-1) 28,26,26
 28   JO=JO+1
      IF (JO-KO) 25,25,2
 29   continue
      return
      END ! of subr. sortr
c
      subroutine SPREADd(R,m,spread) ! Urs Kradolfer
c     DIRICHLET SPREAD FUNCTION
c     measures the goodness of the resolution spread, based on the
c     L2 norm of the difference between the resolution matrix and the
c     identity matrix [ Dirichlet spread function ]
c     R = I   <==>  spread(R) = 0
c
      implicit none
      integer m
      real r(m,m)
      real spread
      integer i,j
c
      spread=0.0
      do i=1,m
         do j=1,m
            if(i.eq.j) spread=spread+(r(i,j)-1)*(r(i,j)-1)
            if(i.ne.j) spread=spread+r(i,j)*r(i,j)
         enddo
      enddo
      return
      end ! of subr. SPREADd
c
      subroutine SPREADb(R,m,spread) ! Urs Kradolfer
c     BACKUS-GILBERT SPREAD FUNCTION
c     measures the goodness of the resolution spread, based on the
c     L2 norm of the WEIGHTED difference between the resolution matrix and the
c     identity matrix [ Backus-Gilbert spread function ]
c     R = I   <==>  spread(R) = 0
c
      implicit none
      integer m
      real r(m,m)
      real spread
      integer j,i
c
      spread=0.0
      do i=1,m
         do j=1,m
            spread=spread+(i-j)*(i-j)*r(i,j)*r(i,j)
         enddo
      enddo
      return
      end ! of subr. spreadb
c
      subroutine SVDSOLUK(Ain,Bin,m,eigmin,X,S,ale,COV,R)
c
c     Urs Kradolfer, 24. April 1987
c
c     solve NORMAL EQUATIONS  Gt*G * X = Gt*R
c                               A  * X  =  B    , so A is symmetric !!!
c
c     Input : left-hand-side A(mxm), right-hand-side  B(m),
c             eigmin : >=0  cutoff-value for solution
c                      = -1 only S (singular values) are returned
c
c     Output: solution-vector X(m), singular-(eigen-)values S(m),
c             diagonal elements of covariance-matrix COV(m),
c             resolution-matrix R(mxm)
c             Average Logarithmic Eigenvalue ALE :
c
c             ALE = 1./(m-izero) * sum(i=1...m) log10(S(i))
c
c             where izero = # zero-singular-values
c
c             before, LGCN was used !!!
c             GCN = LOG10 [ average quotient of 1stEigenvalue to ithEigenvalue
c                           / average eigenvalue ]
c
      implicit none
      integer m
      real eigmin,ale,sum
      integer mm,mm2,i,j,ier,izero,nfre
      parameter (mm=4,mm2=8)
      real Ain(mm,mm),Bin(mm),A(4,4),B(mm),
     &     X(mm),S(mm),COV(mm),R(mm,mm)
      real At(mm,mm)
      real WK(mm2,2)  ! must be dimensioned to : 2m,2
c
      if(m.gt.4)stop'SVDSOLUK>>> m.gt.4 ; redimension arrays !!!'
c
c     copy input-matrix/-vector into local ones (so the input doesn't change) :
c
      do i=1,m
         do j=1,m
            a(i,j)=ain(i,j)
         enddo
         b(i)=bin(i)
         s(i)=0.0
         cov(i)=0.0
      enddo
      do i=1,mm2
         do j=1,2
            wk(i,j)=0.0
         enddo
      enddo
c
c     do SVD (singular value decomposition) of matrix A = V * S * Vt   :
c
      call LSVDF(A,m,m,m,B,m,1,S,WK,ier)   ! IMSL-routine
c
c     determine ALE  ( Average Logarithmic Eigenvalue )
c
c      izero=0
c      if(m.eq.1)then
c         GCN=1.0
c      else
c         GCN=0.0
c         sum=s(1)
c         do i=2,m
c            if(s(i).gt.0.0)then
c               GCN=GCN+s(1)/s(i)
c               sum=sum+s(i)
c            else
c               izero=izero+1
c            endif
c         enddo
c      endif
c      sum=sum/m   ! average eigenvalue.
c      GCN=izero*100.+log10((GCN/float((m-1)-izero))/sum)
cc      ale=0.0
      izero=0
cc      do i=1,m
cc            ale=ale+log10((s(i)/s(1)))
ccc            ale=ale+log10((s(i)/s(1)) + 1e-10)
cc      enddo
cc      ale=-ale/float(m)
      call ALEsubr(s,m,ale)
c
      if(eigmin.lt.0) return
c
c     determine degree of freedom (neglect sing. vals. < cutoff-value eigmin) :
c
      nfre=0
      do i=1,m
         x(i)=0.0
         if(s(i).gt.eigmin)then
            nfre=nfre+1
         else
ccc            write(6,*)'SVDSOLUK>>> skipping singular value = ',s(i)
            s(i)=0.0        ! set eigenvalue to zero
            do j=1,m
               a(j,i)=0.0   ! set appr. COLUMN of V (=eigenvector!!!) to zero
            enddo
         endif
      enddo
c
c     calculate solution-vector X = V * inv(S) * Vt*B    :
c                            [  X = A * inv(S) *   B   ]
c
      do i=1,m
         do j=1,nfre
            x(i)=x(i)+a(i,j)*b(j)/s(j)
         enddo
      enddo
c
c     calculate diagonal elements of UNSCALED covariance matrix :
c
c     C = datvar * V * inv(inv(S)) * Vt
c
c     C(i,j)= SUM(k=1...m)  V(i,k)*V(j,k)/(S(k,k)*S(k,k))  <-- for all elements
c     C(i,i)= SUM(k=1...m)  V(i,k)*V(i,k)/(S(k,k)*S(k,k))  <-- diag. elem. only
c
      do i=1,m
         sum=0.0
         do j=1,nfre
ccc            write(6,*)'i=',i,'j=',j,'V(i,j)=',a(i,j)
            sum=sum+a(i,j)*a(i,j)/(s(j)*s(j))
         enddo
         cov(i)=sum
      enddo
c
c     compute resolution-matrix:  R = V * Vt    :
c                               [ R = A * At ]
c
      call MATRTRAN(A,m,m,At)
      call MATRMULT(A,m,m,At,m,m,R,m,m)
c
      return
      end ! of subr. SVDSOLUK
c
cek     end of vel_math.f
c
cek     begin of vel_ray.f
c
      subroutine RAYPOINTCHECK(rp,nrp,staname)
c
c     checks, whether all the raypoints with z<0 are below the actual
c     surface of the topography (of Switzerland).
c
      implicit none
      include 'vel_com.f'
c
      integer nrp
      real rp(3,inrpmax)
      character*6 staname
c
      integer j
      real zzz,dzzz
      do j=1,nrp
         if(rp(3,j).lt.0.0)then
            call CHTOP(-rp(1,j),rp(2,j),zzz,topo1file,topo2file)
            dzzz=rp(3,j)-zzz
            if(dzzz.lt.0.0)then
               write(6,'(1x,''ray in the air... ! rp3='',f6.3,
     &         '' ZZ='',f6.3,'' dz='',f6.3,'' rp# ='',i2,'' nrp='',i2,
     &         '' STN='',a6,''i '',i4)')
     &         rp(3,j),zzz,dzzz,j,nrp,staname,isingle
            endif
         endif
      enddo
c
      RETURN
      end ! of subr. raypointcheck
c
      subroutine LAYERHIT(rp,nrpdeep,nl,nrp,mll)
c
c     counts the number of hits in each layer.
c
      implicit none
      include 'vel_com.f'
c
      integer nrpdeep,nl,nrp,mll
      real rp(3,inrpmax)
      real rpmax,avraydepth
      integer j,jlay
c
c
c     find deepest layer hit by this ray (= refraction layer or eventlayer
c                                           for direct wave)
      rpmax=-999.
      do j=1,nrp
         if(rp(3,j).gt.rpmax)then
            rpmax=rp(3,j)  ! = max. depth of ray
            nrpdeep=j   ! index of deepest raypoint
         endif
      enddo
      rpmax=rpmax+0.000001   ! avoid raypoint on layer-boundary
      lmax=1
      do j=2,nl
         if(h(j).gt.rpmax)then   !  h(j)= top of layer j [km]
            lmax=j-1  ! lmax = layer number for depth rpmax
            goto 4000
         endif
      enddo
      lmax=nl  ! lowest layer
 4000 continue
      if(rp(3,nrpdeep).eq.rp(3,nrpdeep+1))then   ! ray is horizontal
         irefrlayer(lmax)=irefrlayer(lmax)+1       ! headwave in layer LMAX
         refraylen(lmax)=refraylen(lmax)+
     &   SQRT( (rp(1,nrpdeep)-rp(1,nrpdeep+1))**2 +
     &         (rp(2,nrpdeep)-rp(2,nrpdeep+1))**2 )
      else
         if(mll.eq.0)then
            noheadwave=noheadwave+1
         else
            irefllayer(lmax-1)=irefllayer(lmax-1)+1
c                                           ! -1 because of rpmax=rpmax+0.000001
         endif
      endif
      avhraylen=avhraylen+
     & SQRT( (rp(1,1)-rp(1,nrp))**2 +
     &       (rp(2,1)-rp(2,nrp))**2 )
      avvraylen=avvraylen+
     & ABS( rp(3,nrpdeep)-rp(3,nrp) )
c
      do j=2,nrp
         avraydepth=( rp(3,j)+rp(3,j-1) ) / 2.
         avraydepth=avraydepth+0.000001    ! avoid raydepth on layer-boundary
         do jlay=1,nl
            if(avraydepth.ge.h(jlay).and.
     &         avraydepth.lt.(h(jlay)+thk(jlay)) )then
               hitlay(jlay,1)=hitlay(jlay,1)+1
               hitlay(jlay,2)=hitlay(jlay,2)+
     &                        SQRT( (rp(1,j)-rp(1,j-1))**2 +
     &                              (rp(2,j)-rp(2,j-1))**2 )
               hitlay(jlay,3)=hitlay(jlay,3)+
     &                        ABS( rp(3,j)-rp(3,j-1) )
            endif
         enddo
      enddo
      RETURN
      end ! of subr. layerhit
c
      subroutine CHECKRAYPATH(rp,nrp)
c
c
      implicit none
      include 'vel_com.f'
c
c     do not allow two identical raypoints !
c     (first two raypoints may be identical; e.g. if hypocenter is on
c      layer-boundary)
c
      real rp(3,inrpmax)
      integer nrp
c
      if(rp(1,1).eq.rp(1,2).and.
     &   rp(2,1).eq.rp(2,2).and.
     &   rp(3,1).eq.rp(3,2))then
chrm         do j=1,nrp-1
chrm            rp(1,j)=rp(1,j+1)          ! delete first ray-point
chrm            rp(2,j)=rp(2,j+1)          ! and move elements in array RP 'down'
chrm            rp(3,j)=rp(3,j+1)
chrm         enddo
chrm         nrp=nrp-1                     ! nr_of_raypoints is now smaller by 1 !!
chrm     move hypocenter 10 cm away from the layer boundary (upwards)
         rp(3,1) = rp(3,1) - 0.0001        	                                
         RETURN
      endif
c
      RETURN
      end ! of subr. checkraypath
c
      subroutine BENDRAY(rp,nrp,staname,vtop,ttt)
c
      implicit none
      include 'vel_com.f'
c
      real rp(3,inrpmax)
      integer nrp
      real vtop,ttt
      character*6 staname
c
      logical hypocintop
      real rpn(3,10)
      real xydist,takeoff_angle,arrive_angle
      real deltat1,deltat2,dx1,dy1,dz1,xyz1
      real dx2,dy2,dz2,xyz2,zzz,ttt1new,xyz1n
      real ttt2new,xyz2n,ttt1old,ttt2old
      integer nrpm1,jnrp,jnrpm1
c
      if(icoordsystem.ne.2) RETURN   ! topo-array only for switzerland available
c
      xydist=
     &   SQRT( (rp(1,nrp)-rp(1,1))**2
     &        +(rp(2,nrp)-rp(2,1))**2 )
      if(xydist.gt.10.0) RETURN !epicentral-distance too big for 'airy' rays!(?)
c
      takeoff_angle=
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2 )
     &   /
     &   SQRT( (rp(1,2)-rp(1,1))**2
     &        +(rp(2,2)-rp(2,1))**2
     &        +(rp(3,2)-rp(3,1))**2 )
      takeoff_angle=57.296*ASIN(takeoff_angle)
      if( (rp(3,2)-rp(3,1)) .lt. 0.0 )then
         takeoff_angle=180.-takeoff_angle     ! ray is going upwards
      endif
c
      arrive_angle=
     &   SQRT( (rp(1,nrp)-rp(1,(nrp-1)))**2
     &        +(rp(2,nrp)-rp(2,(nrp-1)))**2 )
     &   /
     &   SQRT( (rp(1,nrp)-rp(1,(nrp-1)))**2
     &        +(rp(2,nrp)-rp(2,(nrp-1)))**2
     &        +(rp(3,nrp)-rp(3,(nrp-1)))**2 )
      arrive_angle=57.296*ASIN(arrive_angle)
      if( (rp(3,2)-rp(3,1)) .lt. 0.0 )then
         arrive_angle=180.-arrive_angle     ! ray is coming upwards
      endif
c
c     write(6,'(1x,''takeoff_angle='',f6.1,''  arrive_angle='',
c    &          f6.1)') takeoff_angle, arrive_angle
c
      deltat1=0.0
      deltat2=0.0
      hypocintop=.false.
      if(rp(3,1).lt.0.0)then
         hypocintop=.true.
c        compute XYZ1: the length between the first 2 raypoints
         dx1=rp(1,2)-rp(1,1)
         dy1=rp(2,2)-rp(2,1)
         dz1=rp(3,2)-rp(3,1)
         xyz1=SQRT(dx1**2 + dy1**2 + dz1**2)
      endif
c
c     if takeoff_angle at surface < 45  --> ray too steep to be 'airy'
c
      if(hypocintop.and.takeoff_angle.le.45.0) hypocintop=.false.
      if(hypocintop)then
         continue   ! stay here; ray-start could be in the air...
      else
         if(arrive_angle.ge.135.0) RETURN ! ray too steep to be 'airy'
      endif
c
      nrpm1=nrp-1
c     compute XYZ2: the length between the last 2 raypoints
      dx2=rp(1,nrp)-rp(1,nrpm1)
      dy2=rp(2,nrp)-rp(2,nrpm1)
      dz2=rp(3,nrp)-rp(3,nrpm1)
      xyz2=SQRT(dx2**2 + dy2**2 + dz2**2)
c
c     average (3D-)ray-length in the top-layer is about 2 kilometers;
c     now split ray in top-layer into 4 parts of about equal length:
      if(hypocintop)then
         dx1=dx1/4.
         dy1=dy1/4.
         dz1=dz1/4.
      endif
      dx2=dx2/4.
      dy2=dy2/4.
      dz2=dz2/4.
      if(hypocintop)then
         do jnrp=1,4  ! first raypoint = first old raypoint!
            rpn(1,jnrp)=rp(1,1)+(jnrp-1)*dx1
            rpn(2,jnrp)=rp(2,1)+(jnrp-1)*dy1
            rpn(3,jnrp)=rp(3,1)+(jnrp-1)*dz1
         enddo
         jnrp=5   ! set last raypoint to second old raypoint to keep accuracy!
         rpn(1,jnrp)=rp(1,2)
         rpn(2,jnrp)=rp(2,2)
         rpn(3,jnrp)=rp(3,2)
      endif
      do jnrp=6,9  ! first raypoint = second-last old raypoint!
         rpn(1,jnrp)=rp(1,nrpm1)+(jnrp-6)*dx2
         rpn(2,jnrp)=rp(2,nrpm1)+(jnrp-6)*dy2
         rpn(3,jnrp)=rp(3,nrpm1)+(jnrp-6)*dz2
      enddo
      jnrp=10 ! set last raypoint to last old raypoint to keep accuracy!
      rpn(1,jnrp)=rp(1,nrp)
      rpn(2,jnrp)=rp(2,nrp)
      rpn(3,jnrp)=rp(3,nrp)
c
c     now set raypoints min. 100meters below surface (if necessary):
c
      if(hypocintop)then
         do jnrp=2,5
            call CHTOP(-rpn(1,jnrp),rpn(2,jnrp),zzz,
     &                 topo1file,topo2file)
            zzz=zzz+0.1  ! 100meters below surface
            if(rpn(3,jnrp).lt.zzz) rpn(3,jnrp)=zzz  ! set hypoc. below surface !
         enddo
      endif
      do jnrp=6,9
         call CHTOP(-rpn(1,jnrp),rpn(2,jnrp),zzz,
     &              topo1file,topo2file)
         zzz=zzz+0.1  ! 100meters below surface
         if(rpn(3,jnrp).lt.zzz) rpn(3,jnrp)=zzz   ! set hypoc. below surface !
      enddo
c
c     calculate ray-length-difference and correct traveltime:
c
      ttt1new=0.0
      if(hypocintop)then
         do jnrp=2,5
            jnrpm1=jnrp-1
            xyz1n=SQRT( (rpn(1,jnrp)-rpn(1,jnrpm1))**2
     &                 +(rpn(2,jnrp)-rpn(2,jnrpm1))**2
     &                 +(rpn(3,jnrp)-rpn(3,jnrpm1))**2 )
            ttt1new=ttt1new+xyz1n/vtop
         enddo
      endif
      ttt2new=0.0
      do jnrp=7,10
         jnrpm1=jnrp-1
         xyz2n=SQRT( (rpn(1,jnrp)-rpn(1,jnrpm1))**2
     &              +(rpn(2,jnrp)-rpn(2,jnrpm1))**2
     &              +(rpn(3,jnrp)-rpn(3,jnrpm1))**2 )
         ttt2new=ttt2new+xyz2n/vtop
      enddo
      if(hypocintop)then
         ttt1old=xyz1/vtop
         deltat1=(ttt1new-ttt1old)
      endif
      ttt2old=xyz2/vtop
      deltat2=(ttt2new-ttt2old)
      if(ABS(deltat1).gt.1e-5.or.ABS(deltat2).gt.1e-5)then
         write(6,'(1x,''BENDRAY>>> ray bended below surface!'',
     &                '' Station: '',a6)') staname
         write(6,'(1x,''dt1='',f6.3,''   dt2='',f6.3)') deltat1,deltat2
      endif
      ttt=ttt+deltat1+deltat2
c
      RETURN
      end ! of subr. bendray
c
      subroutine RAYPATH(nx,ny,nz,x,y,z,vel,nl,thk,d,v,vsq,xe,ye,ze,
     & xr,yr,zr,rp,nrp,nrtn,jl,tkj,itype,ttt,MLL,sterr,direrr,refrerr,
     & reflerr,DTDDrefl,DTDHrefl)
c
c       raypath attempts to approximate the least time path for a
c  seismic ray traveling in a medium with velocities specified at
c  points on a 3-dimensional grid.  for very small delta a straight
c  line path is chosen.  for larger delta a layered velocity model
c  approximating the region between the event and the receiver is
c  constructed and the least time path in that layered model is chosen
c  both direct and refracted rays are considered and the possibility
c  of low velocity layers is allowed.
c       the 3 - dimensional grid should be specified by the
c  intersections of 3 sets of planes normal to the x,y and z coordinate
c  axes.  the intersection of the planes with the axes are specified
c  by x(i), y(j), and z(k), and the spacing of these may be
c  completely arbitrary.
c       subr. layers produces layer velocities from the
c  velocities of the 3 - dimensional model between the event and
c  the receiver, and because these will differ for every event -
c  receiver pair, layers must be called every time raypath is
c  called (unless delta is very small).  the number, thicknesses,
c  and depths of the layers should depend only on the array z(k),
c  and need only to be calculated once for a given velocity model.
c  this is why nl, thk(l), and d(l) are included as input parameters
c  for raypath.  when calling layers, the event and receiver must
c  be located within the horizontal confines of the grid, that is,
c  x(1) <= xe, xr <= x(nx) and y(1) <= ye, yr <= y(ny).
c       subrs. refract and direct determine refracted and direct
c  ray travel times and other parameters in a relatively standard
c  fashion, but are sufficiently general to account for low velocity
c  layers.
c       once the nature of the fastest ray  has been determined, the
c  calculation of a set of raypoints along the path is straightforward
c  and is performed by subr. strpath, refpath, or dirpath.  these
c  points consist of the event and receiver coordinates and, in the cases
c  of refpath and dirpath, every intermediate point at which the ray
c  crosses a boundary between layers.
c       because of the extensive branching in raypath, a parameter
c  to index the the 7 return statements in the main subr.,
c  nrtn is included and takes on values from 1 to 7 (to 8 with extension
c  of reflected waves ! U.K. 11/1986).
c
c  input:
c           nx,ny,nz - number of grid planes in x,y,z directions
c     x(i),y(j),z(k) - gridpoint coordinates
c         vel(i,j,k) - velocity at grid point x(i),y(j),z(k)
c
c                 nl - number of layers (should equal nz)
c             thk(l) - thickness of layer l
c               d(l) - depth to top of layer l
c
c           xe,ye,ze - event coordinates
c           xr,yr,zr - receiver coordinates
c
c  output:
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c               nrp - number of raypoints
c    rp(1,2,or 3,i) - coordinates of raypoint i
c
c  important internal arrays and variables:
c
c            delta - horizontal distance between source and receiver
c            depth - depth of event below receiver
c               jl - number of event layer
c              tkj - depth of event from top of event layer
c
c               kk - refracting layer for fastest refracted ray
c           didjkk - critical distance for refraction in layer kk
c             tref - refracted ray travel time
c           xovmax - an upper bound on delta for which the direct ray
c                       can be the first arrival
c
c           salpha - sine of the takeoff angle of the direct ray
c            deljl - horizontal travel distance of the
c                        direct ray in the event layer
c             tdir - direct ray travel time
c             nrtn - index to point of return from raypath
c
c  subrs. called:
c
c          layers      - constructs  averaged layered model
c
c          stpath      - straight raypath
c          strpath     - NOT USED !! straight raypath from event to receiver
c
c          direct      - for the direct ray:
c                          salpha, deljl, tdir
c          dirpath     - direct ray path
c
c          reflect1     - for reflected ray
c          reflect     - NOT USED !! for reflected ray
c          reflectpath - reflected ray path
c
c          refract     - for the fastest refracted ray:
c                         kk, didjkk, tref, xovmax
c          tiddid      - travel-time intercept & crit. distance for a ray
c          refpath     - refracted ray path
c
      implicit none
      integer nx,ny,nz,nl,nrp,nrtn,jl,itype,mll
      real xe,ye,ze,xr,yr,zr,tkj,ttt,sterr
      real direrr,refrerr,reflerr,dtddrefl,dtdhrefl
      integer l,ierr,kk
      real depth,delta,tdir,salpha,deljl,trefl
      real ain,tref,didjkk,xovmax
      real x(nx),y(ny),z(nz),vel(nx,ny,nz)
      real d(nl),thk(nl),v(nl),vsq(nl)
cek max. nr of raypoints =50
      real rp(3,200)
c
1313  continue   ! restart here in case of an error in REFLECT !!
c
      sterr=0.0
      direrr=0.0
      refrerr=0.0
      reflerr=0.0
c
      if(ze.lt.d(1)) stop'RAYPATH>>> earthquake ABOVE model!'
      if(zr.lt.d(1)) stop'RAYPATH>>> station ABOVE model!'
c
      depth=ze-zr
      delta=sqrt((xr-xe)**2+(yr-ye)**2)
      l=nl
      continue
23107 if(.not.(d(l).ge.ze))goto 23108
c
c   determine event layer
c
      l=l-1
      goto 23107
23108 continue
      jl=l
      d(1)=zr
c   adjust for nonzero depth of reciever
c
      thk(1)=d(2)-zr
      tkj=ze-d(jl)
      if(.not.(delta.lt..05))goto 23105
      call stpath(xe,ye,ze,xr,yr,zr,rp,nrp,ttt,nl,jl,tkj,v,d,thk,sterr)
c
c   if delta is small, set a straight path
c
      nrtn=1
      return
23105 continue
c
c   assign averaged layer velocities
c
      if(itype.eq.1) go to 1
      call layers(nx,ny,nz,x,y,z,vel,xe,ye,xr,yr,nl,v,vsq)
1     continue
      if(.not.(jl.eq.nl))goto 23109
c
c   consider only the direct ray if the event is in the half space
c
      call direct1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,salpha,deljl)
      call dirpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,jl,tkj,
     2 salpha,deljl,rp,nrp,direrr)
      nrtn=2
      ttt=tdir
      return
23109 continue
c
c
      if(mll.le.0)goto 23110
c     MLL > 0  : a reflection layer is specified; ray MUST be reflected!
ctest      call reflect(nl,v,vsq,thk,jl,tkj,delta,depth,mll,trefl,ain,ierr,
ctest     &             DTDDrefl,DTDHrefl)
      call reflect1(nl,v,vsq,thk,jl,tkj,delta,depth,mll,trefl,ain,ierr,
     &             DTDDrefl,DTDHrefl)
      if(ierr.ne.0)then
         write(6,*)'trying another ray-type...'
         write(16,*)'trying another ray-type...'
         MLL=0      ! maybe it works with another ray-type...
         goto 1313  ! retry
      endif
      call reflectpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,jl,tkj,
     &                 ain,mll,rp,nrp,reflerr)
c
      ttt=trefl
      nrtn=8   ! reflected wave
      RETURN
23110 continue
c
c   find the first refracted ray to arrive
c
      call refract(nl,v,vsq,thk,jl,tkj,delta,kk,tref,didjkk,xovmax)
      if(.not.(delta.gt.xovmax))goto 23111
c   for delta greater than xovmax
c   1st arrival is the refracted ray
c
      call refpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,tkj,jl,
     2 kk,didjkk,rp,nrp,refrerr)
      nrtn=3
      ttt=tref
      return
23111 continue
      if(.not.(jl.eq.1))goto 23113
      tdir=(sqrt(tkj**2+delta**2))/v(1)
      if(.not.(tref.lt.tdir))goto 23115
c
c   if refracted ray is 1st arrival
c
      call refpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,tkj,jl,
     2 kk,didjkk,rp,nrp,refrerr)
      nrtn=4
      ttt=tref
      return
23115 continue
c
c   if direct ray is 1st arrival
c
      call stpath(xe,ye,ze,xr,yr,zr,rp,nrp,ttt,nl,jl,tkj,v,d,thk,sterr)
      nrtn=5
      ttt=tdir
      return
23113 continue
c
c   for event below layer 1
c
      call direct1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,salpha,deljl)
      if(.not.(tref.lt.tdir))goto 23117
c
c   if refracted ray is 1st arrival
c
      call refpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,tkj,jl,
     2 kk,didjkk,rp,nrp,refrerr)
      ttt=tref
      nrtn=6
      return
23117 continue
c
c   if direct ray is 1st arrival
c
      call dirpath(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,jl,tkj,
     2 salpha,deljl,rp,nrp,direrr)
      nrtn=7
      ttt=tdir
      return
      end ! of subr. raypath
c
      subroutine LAYERS (nx,ny,nz,x,y,z,vel,xe,ye,xr,yr,nl,v,vsq)
c
c  Subr. layers converts a 3 - dimensional velocity model to a
c  1 - dimensional layered model by averaging slownesses between
c  specified event and receiver positions.  the input
c  model must consist of velocities specified at grid
c  points at the intersections of 3 orthogonal grid planes.
c  slownesses are averaged along a segment joining the
c  horizontal event and receiver positions in each grid plane
c  specified by z(k).
c
c  input:
c           nx,ny,nz - numbers of grid planes
c     x(i),y(j),z(k) - grid point coordinates
c         vel(i,j,k) - velocity at point x(i),y(j),z(k)
c              xe,ye - horizontal event coordinates
c              xr,yr - horizontal receiver coordinates
c                 nl - number of layers (should equal nz)
c
c  output:
c           v(l) - velocity of layer l
c         vsq(l) = v(l) ** 2
c
      implicit none
      integer nx,ny,nz,nl
      real xe,ye,xr,yr
      real x(nx),y(ny),z(nz),vel(nx,ny,nz),v(nl),vsq(nl)
      real w(4),u(100)
      integer i,ie,j,je,ir,jr,imin,jmin,numseg,numpts,l
      integer ipts,ipoint,jpoint
      real a,b,c,d,denom,vpt1,vpt2,vpt3,vpt4,vpoint
      real xseg,yseg,xpoint,ypoint
c
c   locate the event and the reciever within horizontal
c   grid rectangles
c
      i=1
      continue
23119 if(.not.(x(i).lt.xe))goto 23120
      i=i+1
      goto 23119
23120 continue
      ie=i-1
      if(.not.(ie.eq.0))goto 23121
      ie=1
23121 continue
      j=1
      continue
23123 if(.not.(y(j).lt.ye))goto 23124
      j=j+1
      goto 23123
23124 continue
      je=j-1
      if(.not.(je.eq.0))goto 23125
      je=1
23125 continue
      i=1
      continue
23127 if(.not.(x(i).lt.xr))goto 23128
      i=i+1
      goto 23127
23128 continue
      ir=i-1
      if(.not.(ir.eq.0))goto 23129
      ir=1
23129 continue
      j=1
      continue
23131 if(.not.(y(j).lt.yr))goto 23132
      j=j+1
      goto 23131
23132 continue
      jr=j-1
      if(.not.(jr.eq.0))goto 23133
      jr=1
23133 continue
      imin=min0(ie,ir)
      jmin=min0(je,jr)
c
c   choose the number of points to use for averaging slownesses
c   and calculates x and y increments between points
c
      numseg=iabs(ir-ie)+iabs(jr-je)+1
      numpts=numseg+1
      xseg=(xr-xe)/(numseg*1.0)
      yseg=(yr-ye)/(numseg*1.0)
      do 23135l=1,nl
      u(l)=0
23135 continue
23136 continue
c
c   loop to sum slownesses in each layer, u(l), at points
c   between the event and the reciever
c
      do 23137ipts=1,numpts
      xpoint=xe+(ipts-1)*xseg
      ypoint=ye+(ipts-1)*yseg
c
c   locate (xpoint,ypoint) within a grid rectangle
c
      i=imin
      continue
23139 if(.not.(xpoint.gt.x(i)))goto 23140
      i=i+1
      goto 23139
23140 continue
      ipoint=i-1
      if(.not.(ipoint.eq.0))goto 23141
      ipoint=1
23141 continue
      j=jmin
      continue
23143 if(.not.(ypoint.gt.y(j)))goto 23144
      j=j+1
      goto 23143
23144 continue
      jpoint=j-1
      if(.not.(jpoint.eq.0))goto 23145
      jpoint=1
23145 continue
c
c   choosing weighting factors for interpolation within
c   a horizontal rectangle on the velocity grid
c
      a=xpoint-x(ipoint)
      b=x(ipoint+1)-xpoint
      c=ypoint-y(jpoint)
      d=y(jpoint+1)-ypoint
      denom=(a+b)*(c+d)
      w(1)=b*d/denom
      w(2)=a*d/denom
      w(3)=b*c/denom
      w(4)=a*c/denom
c
c   loop to interpolate for the velocity, vpoint, at
c   (xpoint,ypoint) in each layer and add its
c   reciprocal to the sum of the slownesses in that layer
c
      do 23147l=1,nl
      vpt1=w(1)*vel(ipoint,jpoint,l)
      vpt2=w(2)*vel(ipoint+1,jpoint,l)
      vpt3=w(3)*vel(ipoint,jpoint+1,l)
      vpt4=w(4)*vel(ipoint+1,jpoint+1,l)
      vpoint=vpt1+vpt2+vpt3+vpt4
      u(l)=u(l)+1/vpoint
23147 continue
23148 continue
23137 continue
23138 continue
c
c   loop to calculate the velocity from the sum of the
c   slownesses for each layer
c
      do 23149l=1,nl
      v(l)=(numpts*1.0)/u(l)
      vsq(l)=v(l)**2
23149 continue
23150 continue
      return
      end ! of subr. layers
c
      subroutine STPATH(xe,ye,ze,xr,yr,zr,rp,nrp,tt,nl,jl,tkj,v,d,thk,
     &                  sterr)
      implicit none
      real xe,ye,ze,xr,yr,zr,tt,tkj,sterr
      integer nrp,nl,jl
c
      integer jl1,j,j1,jb,jt
      real v(nl),d(nl),thk(nl)
c
      real rp(3,200)
      nrp=jl+1
c  assign points along raypath
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
c
      if (nrp.eq.2) go to 50
c
      jl1=jl-1
cuk      do 25 j=1,jl1    ! calculate up to the top and then compute an error
      do 25 j=1,jl1+1
      j1=j+1
      rp(1,j1)=xe
      rp(2,j1)=ye
      jb=jl-j+1
      rp(3,j1)=d(jb)
   25 continue
c
   50 continue
      if(nrp.gt.2)then
         sterr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &               (rp(3,nrp)-zr)**2 ) *1000.
      endif
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
c
c  compute travel time
      tt=tkj/v(jl)
      if (nrp.eq.2) return
c
      do 75 j=1,jl1
      jt=jl-j
      tt=tt+thk(jt)/v(jt)
   75 continue
      return
      end ! of subr. stpath
c
      subroutine STRPATH (xe,ye,ze,xr,yr,zr,rp,nrp)
c
c       strpath assigns raypoint coordinates for a straight path
c  between the event and the receiver.
c
c  input:  (xe,ye,ze) - event coordinates
c          (xr,yr,zr) - receiver coordinates
c
c  output:   rp(1,2,or 3,i) - coordinates of raypoint i
c                       nrp - number of raypoints = 2
c
      implicit none
      real xe,ye,ze,xr,yr,zr
      integer nrp
      real rp(3,2)
      nrp=2
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
      rp(1,2)=xr
      rp(2,2)=yr
      rp(3,2)=zr
      return
      end ! of subr. strpath
c
      subroutine DIRECT1(nl,v,vsq,thk,jl,tkj,delta,depth,tdir,u,x)
c
c       for the direct seismic ray from an event to a receiver in
c  a layered velocity structure, direct predicts the travel time, the
c  sine of the takeoff angle, and the horizontal distance of travel in
c  the event layer.  the receiver must be located at the top of layer
c  1 and the event must be located below layer 1.  low velocity
c  layers are permitted.
c       to find the takeoff angle of the ray, a numerical approach
c  is required.  the basic scheme adopted here is the method of false
c  position.  (see acton, 1970, 'numerical methods that work,' for
c  example.)  first, the characteristics of the fastest layer
c  between the event and the surface are determined.  these permit
c  placing definite lower and upper bounds, ua and ub, on the
c  sine of the takeoff angle.  in turn, lower and upper bounds, xa
c  and xb, on the horizontal travel distance in the event layer are
c  determined.  the total horizontal travel distance for a ray with
c  with horizontal travel distance x in the event layer is denoted
c  by del, and the zero of del - delta is obtained by using xa and
c  xb as initial guesses for x in the method of false position
c  from x and tkj, the depth of the event below the top of the event
c  layer, the sine of the takeoff angle, u , is calculated.
c       from u and x, tdir is found by summing the travel time in
c  each layer.  finally, a slight correction to tdir is made, based
c  on the misfit between the final del and delta.
c
c  input:     nl - number of layers
c           v(l) - velocity of layer l
c         vsq(l) = v(l) ** 2
c         thk(l) - thickness of layer l
c             jl - event layer
c            tkj - depth of event from top of event layer
c          delta - horizontal distance between event and receiver
c          depth - depth of event
c
c  output:  tdir - direct ray travel time
c              u - sine of the takeoff angle
c              x - horizontal travel distance in the event layer
c
c
c  find the fastest layer, lmax, above and including jl
c
      implicit none
      integer nl,jl
      real tkj,delta,depth,tdir,u,x
      real v(nl),vsq(nl),thk(nl)
c
      real del,tklmax,vlmax,ua,ub,uasq,ubsq,usq
      real xa,xb,dela,delb,ubdiv,xtest
      integer lmax,j1,l,kount
      del=0.0   ! U.K. 28. Jan. 1987
c
      lmax=jl
      tklmax=tkj
      vlmax=v(jl)
      j1=jl-1
      do 23184l=1,j1
      if(.not.(v(l).gt.vlmax))goto 23186
      lmax=l
      tklmax=thk(l)
      vlmax=v(l)
23186 continue
23184 continue
23185 continue
C CHANGE BY E.KISSLING MARCH 1984
      IF(tklmax.le.0.05) tklmax=0.05
c
c   find initial bounds on the sine of the takeoff angle
c
      ua=(v(jl)/vlmax)*delta/sqrt(delta**2+depth**2)
      ub=(v(jl)/vlmax)*delta/sqrt(delta**2+tklmax**2)
c
c   calculate horizontal travel distances
c
      uasq=ua**2
      ubsq=ub**2
C CHANGE BY E.KISSLING MARCH 1984
      IF(UBSQ.GE.1.) UBSQ=0.99999
      IF(UASQ.GE.1.) UASQ=0.99999
      xa=tkj*ua/sqrt(1.0-uasq)
      if(.not.(lmax.eq.jl))goto 23188
      xb=delta
      goto 23189
23188 continue
      xb=tkj*ub/sqrt(1.0-ubsq)
23189 continue
      dela=xa
      delb=xb
      do 23190l=1,j1
      dela=dela+thk(l)*ua/sqrt(vsq(jl)/vsq(l)-uasq)
      ubdiv=sqrt(vsq(jl)/vsq(l)-ubsq)
      if(ubdiv.GT.1.e-20) GOTO 1002
      write(16,*)'WARNING:'
      write(16,1000) jl,l,lmax,vsq(jl),vsq(l),ubsq,delta,tklmax,vlmax
 1000 format(/,2x,'subr. direct1: ',3i4,2f10.4,f15.6,3f10.5,/)
      ubdiv=1.e-20
 1002 continue
      delb=delb+thk(l)*ub/sqrt(vsq(jl)/vsq(l)-ubsq)
23190 continue
23191 continue    !  NOT used... !!!
c
c   loop to find the zero of del-delta by teh method of false position
c
      do 23192kount=1,25
      if(.not.((delb-dela).lt.0.02))goto 23194
      x=0.5*(xa+xb)
      u=x/sqrt(x**2+tkj**2)
      usq=u**2
      goto 23193
23194 continue
      x=xa+(delta-dela)*(xb-xa)/(delb-dela)
      u=x/sqrt(x**2+tkj**2)
      usq=u**2
      del=x
      do 23196l=1,j1
      del=del+thk(l)*u/sqrt(vsq(jl)/vsq(l)-usq)
23196 continue
23197 continue
      xtest=del-delta
      if(.not.(abs(xtest).lt.0.02))goto 23198
      goto 23193
23198 continue
      if(.not.(xtest.lt.0.0))goto 23200
      xa=x
      dela=del
      goto 23201
23200 continue
      xb=x
      delb=del
23201 continue
23192 continue
23193 continue
c
c   calculate direct ray travel time
c
c
      if(del.eq.0.0) del=x   ! U.K. 28. Jan. 1987
c
      tdir=(sqrt(x**2+tkj**2))/v(jl)
      do 23202l=1,j1
      tdir=tdir+thk(l)*v(jl)/(vsq(l)*sqrt(vsq(jl)/vsq(l)-usq))
23202 continue
23203 continue
      tdir=tdir-(u/v(jl))*(del-delta)
      return
      end ! of subr. direct1
c
      subroutine DIRPATH (xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,
     &jl,tkj,salpha,deljl,rp,nrp,direrr)
c
c       dirpath assigns raypoints with coordinates rp(1,i),rp(2,i),
c  rp(3,i) along a direct seismic raypath in a flat layered earth.
c
c  input:  xe,ye,ze - event coordinates
c          xr,yr,zr - receiver coordinates
c             delta - horizontal distance between event and receiver
c                nl - total number of layers in model
c                     (used only to dimension v vsq and thk)
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c            thk(l) - thickness of layer l
c                jl - event layer
c               tkj - depth of event from top of event layer
c            salpha - sine of takeoff angle
c             deljl - horizontal travel distance in event layer
c
c  output:  rp(1,2,or 3,i) - raypoint coordinates
c                      nrp - number of raypoints
c
      implicit none
      real xe,ye,ze,xr,yr,zr,delta,tkj,salpha,deljl
      real direrr
      integer nl,jl,nrp,nrp1,i,m
      real d1,d2,tng
cek max. nr of raypoints =50
      real v(nl),thk(nl),vsq(nl),rp(3,200)
      d1=(xr-xe)/delta
      d2=(yr-ye)/delta
      nrp=jl+1
c
c   assign raypoint coordinates
c
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
      rp(1,2)=rp(1,1)+deljl*d1
      rp(2,2)=rp(2,1)+deljl*d2
      rp(3,2)=rp(3,1)-tkj
cuk      nrp1=nrp-1
      nrp1=nrp  ! calculate up to the top and then compute an error !
      do 23208i=3,nrp1
      m=nrp-(i-1)
c
c   m is the number of the layer below ray point i
c
      tng=v(m)*salpha/sqrt(vsq(jl)-vsq(m)*salpha**2)
c
c   tng is the tangent of the angle of incidence in layer m
c
      rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
      rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
      rp(3,i)=rp(3,i-1)-thk(m)
23208 continue
23209 continue
      direrr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &             (rp(3,nrp)-zr)**2 ) *1000.
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
      return
      end ! of subr. dirpath
c
      subroutine REFLECT1(nl,v,vsq,thk,jl,tkj,delta,z,mll,trefl,ain,
     &                    ierr,DTDDrefl,DTDHrefl)
c
c     Urs Kradolfer, Nov. 1986
c
      implicit none
c      
      integer nl,ierr
      real tkj,trefl,dtddrefl,dtdhrefl
c     
      integer m,lmax,j1,l,kount
      real del, fac,depth,tklmax,vlmax,ua,ub,uasq,ubsq
      real xa,xb,dela,delb,ubdiv,x,u,usq
      real xtest,tdir
      integer jl   ! hypocenter-layer
      integer mll  ! reflection at bottom of layer MLL
      real thk(nl),v(nl),vsq(nl)
      real div(100)  ! dimension = max. number of layers allowed
      real ain     ! auftauchwinkel in rad
      real z       ! hypocenter-depth
      real delta  ! distance epicenter-receiver
      integer i
c
      ierr=0
      DTDDrefl=0.0
      DTDHrefl=0.0
c
c
      del=0.0   ! U.K. 28. Jan. 1987
c
      if(jl.gt.mll)then   !  hypocenter below reflector !!!
         write(6,*)'hypocenter-layer jl = ',jl
         write(6,*)'reflection at bottom of layer mll = ',mll
         ierr=-1
         write(6,*)'WARNING:'
         write(6,*)'subr. REFLECT1 >>> hypocenter is below reflector!'
         write(16,*)'WARNING:'
         write(16,*)'subr. REFLECT1 >>> hypocenter is below reflector!'
ccc         stop'subr. REFLECT1 >>> hypocenter below reflector!'
      endif
      if(mll.gt.nl)then
         stop'subr. REFLECT1 >>> reflector-nr > number of layers !'
      endif
c
c  determine out of real model (THK) transformed model (DIV)
c
      DO M=1,MLL
C LAYER M IS PASSED ONLY ONCE:
         FAC= 1.
         IF(JL.GT.M) GOTO 505
C LAYER M IS PASSED TWICE
         FAC= 2.
         IF(JL.LT.M) GOTO 505
C LAYER M CONTAINS HYPOCENTER:
         FAC= (2.*THK(M)-tkj)/THK(M)
  505    DIV(M)= FAC*THK(M)
      enddo
      depth=0.0
      do i=1,mll
         depth=depth+div(i)
      enddo
c
c
c  now 'new' model is established; layer-thicknesses are stored in DIV(1...mll)
c
      lmax=mll
      tklmax=div(mll)
      vlmax=v(mll)
      j1=mll-1
      do 23184l=1,j1
      if(.not.(v(l).gt.vlmax))goto 23186
      lmax=l
      tklmax=div(l)
      vlmax=v(l)
23186 continue
23184 continue
23185 continue
C CHANGE BY E.KISSLING MARCH 1984
      IF(tklmax.le.0.05) tklmax=0.05
c
c   find initial bounds on the sine of the takeoff angle
c
      ua=(v(mll)/vlmax)*delta/sqrt(delta**2+depth**2)
      ub=(v(mll)/vlmax)*delta/sqrt(delta**2+tklmax**2)
c
c   calculate horizontal travel distances
c
      uasq=ua**2
      ubsq=ub**2
C CHANGE BY E.KISSLING MARCH 1984
      IF(UBSQ.GE.1.) UBSQ=0.99999
      IF(UASQ.GE.1.) UASQ=0.99999
      xa=div(mll)*ua/sqrt(1.0-uasq)
      if(.not.(lmax.eq.mll))goto 23188
      xb=delta
      goto 23189
23188 continue
      xb=div(mll)*ub/sqrt(1.0-ubsq)
23189 continue
      dela=xa
      delb=xb
      do 23190l=1,j1
      dela=dela+div(l)*ua/sqrt(vsq(mll)/vsq(l)-uasq)
      ubdiv=sqrt(vsq(mll)/vsq(l)-ubsq)
      if(ubdiv.GT.1.e-20) GOTO 1002
      write(16,*)'WARNING:'
      write(16,1000) mll,l,lmax,vsq(mll),vsq(l),ubsq,delta,tklmax,vlmax
 1000 format(/,2x,'subr. reflect1: ',3i4,2f10.4,f15.6,3f10.5,/)
      ubdiv=1.e-20
 1002 continue
      delb=delb+div(l)*ub/sqrt(vsq(mll)/vsq(l)-ubsq)
23190 continue
23191 continue    !  NOT used... !!!
c
c   loop to find the zero of del-delta by teh method of false position
c
      do 23192kount=1,25
      if(.not.((delb-dela).lt.0.02))goto 23194
      x=0.5*(xa+xb)
      u=x/sqrt(x**2+div(mll)**2)
      usq=u**2
      goto 23193
23194 continue
      x=xa+(delta-dela)*(xb-xa)/(delb-dela)
      u=x/sqrt(x**2+div(mll)**2)
      usq=u**2
      del=x
      do 23196l=1,j1
      del=del+div(l)*u/sqrt(vsq(mll)/vsq(l)-usq)
23196 continue
23197 continue
      xtest=del-delta
      if(.not.(abs(xtest).lt.0.02))goto 23198
      goto 23193
23198 continue
      if(.not.(xtest.lt.0.0))goto 23200
      xa=x
      dela=del
      goto 23201
23200 continue
      xb=x
      delb=del
23201 continue
23192 continue
23193 continue
c
c   calculate direct ray travel time
c
c
      if(del.eq.0.0) del=x   ! U.K. 28. Jan. 1987
c
      tdir=(sqrt(x**2+div(mll)**2))/v(mll)
      do 23202l=1,j1
      tdir=tdir+div(l)*v(mll)/(vsq(l)*sqrt(vsq(mll)/vsq(l)-usq))
23202 continue
23203 continue
      tdir=tdir-(u/v(mll))*(del-delta)
c
      trefl=tdir
c---- u is sine of the 'takeoff-angle' in
c---- transformed model = emerging angle at reflector!
c---- Find now takeoff-angle at source!
      ain=u
      if(mll.gt.jl)then
         do i=mll-1,jl,-1
            ain=ain * v(i)/v(i+1)
         enddo
      endif
c
      return
      end ! of subr. reflect1
c
      subroutine REFLECT(nl,v,vsq,thk,jl,tkj,delta,z,mll,trefl,ain,ierr,
     &             DTDDrefl,DTDHrefl)
c
c     Urs Kradolfer, Nov. 1986
c
      implicit none
      integer nl,ierr
      real tkj,trefl,ain,dtddrefl,dtdhrefl
c
      integer nit,m
      real zbig,fac,a,sina,dellit,dddas,sqr,ddda,da,da2
      real tt,dtda,ti,dti,dtdd,sina2,dtdh
      integer jl   ! hypocenter-layer
      integer mll  ! reflection at bottom of layer MLL
      real thk(nl),v(nl),vsq(nl)
      real vsqu(100),vi(100),div(100)  ! dimension = max. number of layers allowed
      real anin    ! auftauchwinkel in rad
      real test4   ! iterate until hypoc-adjustments are less than (TEST4)**-2
      real z       ! hypocenter-depth
      real delta  ! distance epicenter-receiver
      integer i
      real sum
c
      ierr=0
      DTDDrefl=0.0
      DTDHrefl=0.0
c
      test4=0.01                  !  see HYMD01.DAT
      do i=1,nl
         vi(i)=1./v(i)
      enddo
c
c---- determine, which layer contains hypocenter
   10 sum=0.0
      do jl=1,nl
         sum=sum+thk(jl)
         if(z.le.sum)goto 19
      enddo
c---- hypocenter is in the halfspace:
      jl=nl
      tkj=z-sum
      if(jl.eq.0) goto 10
      goto 20
   19 tkj=z-sum+thk(jl)
      goto 20
   20 continue
c     hypocenter is TKJ [km] below top of hypoc-layer
c     JL is the hypocenter layer
c
      if(jl.gt.mll)then   !  hypocenter below reflector !!!
         ierr=-1
         write(6,*)'WARNING:'
         write(6,*)'subr. REFLECT >>> hypocenter is below reflector!'
         write(16,*)'WARNING:'
         write(16,*)'subr. REFLECT >>> hypocenter is below reflector!'
ccc         stop'subr. REFLECT >>> hypocenter below reflector!'
      endif
      if(mll.gt.nl)then
         stop'subr. REFLECT >>> reflector-nr > number of layers !'
      endif
c
C
C  CALCULATE REFLECTED WAVES
C
      NIT= 0
      ZBIG= TKJ
      SUM= 0.0
      DO M=1,MLL
C LAYER M IS PASSED ONLY ONCE:
         FAC= 1.
         IF(JL.GT.M) GOTO 505
C LAYER M IS PASSED TWICE
         FAC= 2.
         IF(JL.LT.M) GOTO 505
C LAYER M CONTAINS HYPOCENTER:
         FAC= (2.*THK(M)-ZBIG)/THK(M)
  505    DIV(M)= FAC*THK(M)
         VSQU(M)= VSQ(MLL)/VSQ(M)
         SUM= SUM+DIV(M)
      enddo
C ITERATE TO FIND EMERGING ANGLE FOR THE GIVEN DISTANCE:
      A= ATAN2(delta,SUM)
  515 SINA= SIN(A)
      if(sina.eq.1.00) sina=0.9999999
      SINA2= SINA*SINA
      DELLIT= 0.
      DDDAS= 0.
      DO M=1,MLL
         SQR= VSQU(M)-SINA2
c
       if(sqr.le.0)then
         write(16,*)'WARNING:                  subr. REFLECT'
         write(16,*)'probably a low velocity layer detected!'
         write(16,*)'sina=',sina
         write(16,*)'sina2=',sina2
         write(16,*)'VSQU(M):  vsqu(',m,')=',vsqu(m)
         write(16,*)'MLL=',mll
         write(16,*)'JL=',jl
         write(6,*)'sina=',sina
         write(6,*)'sina2=',sina2
         write(6,*)'VSQU(M):  vsqu(',m,')=',vsqu(m)
         write(6,*)'MLL=',mll
         write(6,*)'JL=',jl
         write(6,*)
       endif
c
c!!!!!         if(sqr.le.0.) goto 530   !auftauchwinkel gefunden
c   sina2 is always < 1 ;
c   vsqu is always > 1 if no velocity-layers are above the event-layer !!
         DDDA= DIV(M)/SQRT(SQR)
         DELLIT= DELLIT+DDDA
         DDDAS= DDDAS+VSQU(M)*DDDA/SQR
      enddo
      DELLIT= DELLIT*SINA
      DDDAS= DDDAS*COS(A)
      DA= DELLIT-delta
      DA2= DA*DA
      IF(DA2.LE.TEST4) GOTO 530
      NIT= NIT+1
cuk      IF(NIT.LT.15) GOTO 525
      IF(NIT.LT.50) GOTO 525
c      WRITE(6,524) delta
c  524 FORMAT(' REFLECTED WAVE DID NOT CONVERGE WITHIN',
c     1' 50 ITERATIONS AT DISTANCE:',F8.1,
c     2' SET WEIGHT TO 0')
      GOTO 27
  525 A= A-DA/DDDAS
      if(a.ge.1.570796) a=1.5  ! 1.5 ~ 86 deg ;   a = never .gt.  pi/2.  !!!
      GOTO 515
c---- successful iteration:
  530 TT=0.
      DTDA= 0.
      DO M=1,MLL
         TI= DIV(M)*VI(M)/SQRT(1.-SINA2/VSQU(M))
         DTI= TI/(VSQU(M)-SINA2)
         TT= TT+TI
         DTDA= DTDA+DTI
      enddo
      DTDA= DTDA*SINA*COS(A)
      DTDD=  DTDA/DDDAS
      DTDH= -(1.-V(JL)*SINA*DTDD)/(V(JL)*SQRT(1.-SINA2))
      ANIN= SINA
      GOTO 260
C------- IF NOT POSSIBLE, SET P + S WEIGHTS TO ZERO
   27 continue
      write(16,*)'WARNING:'
      write(16,*)'Reflected wave did not converge within 50 iterations!'
      write(16,*)'---> trying another ray-type ...'
      write(6,*)'Reflected wave did not converge within 50 iterations!'
      write(6,*)'sina=',sina,'  delta=',delta
      write(6,*)'da=',da,'  da2=',da2,'  test4=',test4
      write(6,*)'---> trying another ray-type ...'
      write(6,*)
      ierr=50
      return
cc      stop'STOP in subr. REFLECT >>> nit > 50'
c
  260 continue
c---- anin is emerging angle at reflector; find now takeoff-angle at source!
      ain=anin
      if(mll.gt.jl)then
         do i=mll-1,jl,-1
            ain=ain * v(i)/v(i+1)
         enddo
      endif
      trefl=tt
      DTDDrefl=dtdd
      dtdhrefl=dtdh
cc     wain=57.29578*asin(ain)
cc      write(6,*)'takeoff angle = ',wain
      return
c
      end ! of subr. reflect
c
      subroutine REFLECTPATH(xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,
     &jl,tkj,ain,mll,rp,nrp,reflerr)
c
c     Urs Kradolfer, Nov. 1986
c
c       reflectpath assigns raypoints with coordinates rp(1,i),rp(2,i),
c  rp(3,i) along a reflected seismic raypath in a flat layered earth.
c
c  input:  xe,ye,ze - event coordinates
c          xr,yr,zr - receiver coordinates
c             delta - horizontal distance between event and receiver
c                nl - total number of layers in model
c                     (used only to dimension v vsq and thk)
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c            thk(l) - thickness of layer l
c                jl - event layer
c               tkj - depth of event from top of event layer
c               ain - sine of takeoff angle (with respect to downward vertical)
c               mll - reflector is bottom of layer MLL
c
c  output:  rp(1,2,or 3,i) - raypoint coordinates
c                      nrp - number of raypoints
c
c max. nr of raypoints =50
c
      implicit none
c
      real xe,ye,ze,xr,yr,zr,delta,tkj,ain,reflerr
      integer nl,jl,mll,nrp
c
      real d1,d2,salpha,tng,deljl
      integer m,ld,j,i
      real v(nl),thk(nl),vsq(nl),rp(3,200)
c
      d1=(xr-xe)/delta
      d2=(yr-ye)/delta
      nrp=jl+1+(mll-jl)*2+1
c
      salpha=ain
c
c   assign raypoint coordinates
c
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
c      tng=tan(asin(salpha))  ! compiler produces warning with range checking
      tng=asin(salpha)
      tng=tan(tng)       ! must be defined here; in case ld is zero
      deljl=tng*(thk(jl)-tkj) !horiz. travel distance in ev.layer
      rp(1,2)=rp(1,1)+deljl*d1
      rp(2,2)=rp(2,1)+deljl*d2
      rp(3,2)=rp(3,1)+(thk(jl)-tkj)
c---- now arrived at bottom of event-layer
      m=jl
      ld=mll-jl    ! nr of layer to go down to reflector
      j=2
      do i=1,ld
         j=j+1
         m=jl+i
         tng=v(m)*salpha/sqrt(vsq(jl)-vsq(m)*salpha**2)
         rp(1,j)=rp(1,j-1)+thk(m)*tng*d1
         rp(2,j)=rp(2,j-1)+thk(m)*tng*d2
         rp(3,j)=rp(3,j-1)+thk(m)       ! ray is going downward
      enddo
c---- now arrived down at reflector; going upward from now on
      if(m.ne.mll)stop'subr. REFLECTPATH >>> Error in ray-path...'
      j=j+1
      deljl=tng*thk(m)   ! horiz. travel-distance in lowest layer
      rp(1,j)=rp(1,j-1)+deljl*d1
      rp(2,j)=rp(2,j-1)+deljl*d2
      rp(3,j)=rp(3,j-1)-thk(m)
c---- now arrived at top of lowest layer
c---- going up to top of model:
c     do i=j+1,nrp-1   if one would raytrace only to the second-last point
c in order to calculate an error, we calculate up to the receiver:
      do i=j+1,nrp
         m=m-1
c
c   m is the number of the layer below ray point i
c
         tng=v(m)*salpha/sqrt(vsq(jl)-vsq(m)*salpha**2)
c
c   tng is the tangent of the angle of incidence in layer m
c
         rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
         rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
         rp(3,i)=rp(3,i-1)-thk(m)
      enddo
c
      reflerr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &              (rp(3,nrp)-zr)**2 ) *1000.
cc      write(6,*)'REFLECTPATH: raytracer-error was : ',reflerr,' meter'
c---- now assign the correct receiver-coordinates to the last raypoint:
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
      return
      end ! of subr. reflectpath
c
      subroutine REFRACT (nl,v,vsq,thk,jl,tkj,delta,
     &kk,tref,didjkk,xovmax)
c
c       for refracted rays in a layered earth model, refract
c  determines the fastest travel time, tref, the layer
c  in which the fastest ray is refracted, kk, the
c  critical distance for refraction in that layer,
c  didjkk, and an upper bound on delta for which a
c  direct ray can be a first arrival, xovmax.  refract
c  allows for the possibility of low velocity layers.
c       note that there may not be a refracted ray, either because
c  all layers below the event layer are low velocity layers or
c  because for all layers below the event layer which are not low
c  velocity layers the critical distance exceeds delta.  in such
c  cases tref, didjkk, and xovmax are set very large, kk is set to
c  zero, and refract returns to the calling program.
c
c  input:  nl - number of layers
c        v(l) - velocity of layer l
c      vsq(l) - v(l) ** 2
c      thk(l) - thickness of layer l
c          jl - event layer
c         tkj - depth of event in event layer
c       delta - horizontal distance between event and receiver
c
c  output:   kk - refracting layer for fastest refracted ray
c          tref - travel time of fastest refracted ray
c        didjkk - critical distance for refraction in layer kk
c        xovmax - an upper bound on delta for which the direct ray can
c                       be the first arrival
c  internal arrays:
c
c       tr(m) - travel time for refraction in layer m
c     tinj(m) - traveltime intercept
c      tid(m) - terms in travel time intercept which are
c                     independent of tkj
c     didj(m) - critical distance
c      did(m) - terms in critical distance which are
c                     independent of tkj
c
c
c  call subr. tiddid to evaluate tid(m) and
c  did(m), the terms in the travel time intercept and
c  critical distance for a ray refracted in layer m
c  which are independent of tkj.
c
c
c  determine tref, kk, didjkk
c
      implicit none
c
      integer nl,jl,kk,j1,m,lx,m1,l,jx
      real tkj,delta,tref,didjkk,xovmax
      real sqt,tim
      real v(nl),vsq(nl),thk(nl)
      real tid(100),did(100),tinj(100),didj(100),tr(100)
      call tiddid(jl,nl,v,vsq,thk,tid,did)
      tref=100000.
      j1=jl+1
      do 23151m=j1,nl
      if(.not.(tid(m).eq.100000.))goto 23153
      tr(m)=100000.
      goto 23154
23153 continue
      sqt=sqrt(vsq(m)-vsq(jl))
      tinj(m)=tid(m)-tkj*sqt/(v(m)*v(jl))
      didj(m)=did(m)-tkj*v(jl)/sqt
      tr(m)=tinj(m)+delta/v(m)
      if(.not.(didj(m).gt.delta))goto 23155
      tr(m)=100000.
23155 continue
23154 continue
      if(.not.(tr(m).lt.tref))goto 23157
      tref=tr(m)
      kk=m
      didjkk=didj(kk)
23157 continue
23151 continue
23152 continue
c
c   if there is no refracted ray:
c
      if(.not.(tref.eq.100000.))goto 23159
      didjkk=100000.
      xovmax=100000.
      kk=0
      return
23159 continue
c
c   if there is a refracted ray, determine xovmax:
c   find lx, the 1st layer below the event layer which
c   is not a low velocity layer
c
      m=jl+1
      continue
23161 if(.not.(tid(m).eq.100000.))goto 23162
      m=m+1
      goto 23161
23162 continue
      lx=m
c
c   check whether the event is in the 1st layer
c
      if(.not.(jl.eq.1))goto 23163
      xovmax=tinj(lx)*v(lx)*v(1)/(v(lx)-v(1))
      return
23163 continue
      m=jl
c
c   find jx, the 1st layer above and including the event
c   layer which is not a low velocity layer
c
      continue
23165 continue
      tid(m)=0.
      m1=m-1
      do 23168l=1,m1
      if(.not.(vsq(m).le.vsq(l)))goto 23170
      tid(m)=100000.
      goto 23171
23170 continue
      sqt=sqrt(vsq(m)-vsq(l))
      tim=thk(l)*sqt/(v(l)*v(m))
      tid(m)=tid(m)+tim
23171 continue
23168 continue
23169 continue
      m=m-1
c
c  decide whether or not jx=1 and calculate xovmax
c
23166 if(.not.(tid(m+1).lt.100000..or.m.eq.1))goto 23165
23167 continue
      if(.not.(tid(m+1).lt.100000.))goto 23172
      jx=m+1
      xovmax=(tinj(lx)-tid(jx))*v(lx)*v(jx)/(v(lx)-v(jx))
      goto 23173
23172 continue
c
c   jx=1
c
      xovmax=tinj(lx)*v(lx)*v(1)/(v(lx)-v(1))
23173 continue
      return
      end ! of subr. refract
c
      subroutine TIDDID (jl,nl,v,vsq,thk,tid,did)
c
c       tiddid determines the travel time intercept and critical
c  distance for a seismic ray in a layered earth model
c  originating at the top of layer jl, refracting in
c  layer m, and terminating at the top of layer 1.
c
c  input:       jl - event layer
c               nl - number of layers
c             v(l) - velocity of layer l
c           vsq(l) - velocity squared
c           thk(l) - thickness of layer l
c  output:
c           tid(m) - travel time intercept for
c                      refraction in layer m
c           did(m) - critical distance
c
      implicit none
c
      integer jl,nl
c
      integer j1,m,m1,l
      real tid1,tid2,did1,did2,sqt,tim,dimm
      real v(nl),vsq(nl),thk(nl)
      real tid(100),did(100)
      j1=jl+1
      do 23174m=j1,nl
      tid(m)=0.
      did(m)=0.
      tid1=0.
      tid2=0.
      did1=0.
      did2=0.
      m1=m-1
      do 23176l=1,m1
      if(.not.(vsq(m).le.vsq(l)))goto 23178
c
c   if m is a low velocity layer, set tid and did to
c   very large values
c
      tid(m)=100000.
      did(m)=100000.
      goto 23179
23178 continue
      sqt=sqrt(vsq(m)-vsq(l))
      tim=thk(l)*sqt/(v(l)*v(m))
      dimm=thk(l)*v(l)/sqt
      if(.not.(l.lt.jl))goto 23180
c
c   sum for layers above event layer
c
      tid1=tid1+tim
      did1=did1+dimm
      goto 23181
23180 continue
c
c   sum for layers below and including the event layer
c
      tid2=tid2+tim
      did2=did2+dimm
23181 continue
23179 continue
23176 continue
23177 continue
      if(.not.(tid(m).ne.100000.))goto 23182
c
c   calculate tid and did if m is not a low velocity layer
c
      tid(m)=tid1+2*tid2
      did(m)=did1+2*did2
23182 continue
23174 continue
23175 continue
      return
      end ! of subr. tiddid
c
      subroutine REFPATH (xe,ye,ze,xr,yr,zr,delta,nl,v,vsq,thk,
     &tkj,jl,kk,didjkk,rp,nrp,refrerr)
c
c       refpath assigns raypaths with coordinates rp(1,i), rp(2,i),
c  rp(3,i) along a refracted seismic raypath in a flat layered earth.
c
c  input:  xe,ye,ze - event coordinates
c          xr,yr,zr - receiver coordinates
c             delta - horizontal distance between the event
c                          and the receiver
c                nl - total number of layers in model
c                     (used only to dimension v, vsq, and thk)
c              v(l) - velocity of layer l
c            vsq(l) = v(l) ** 2
c            thk(l) - thickness of layer l
c                jl - event layer
c               tkj - depth of event from top of event layer
c                kk - refracting layer
c            didjkk - critical distance for the refracted ray
c
c  output:  rp(1,2,or 3,i) - raypoint coordinates
c                      nrp - number of raypoints
c
cek max. nr of raypoints =50
c
      implicit none
c
      real xe,ye,ze,xr,yr,zr,delta,tkj,didjkk,refrerr
      integer nl,jl,kk,nrp
c
      real thkjl,d1,d2,tng
      integer nrpd,i,m,nrp1,nrpd2
      real v(nl),thk(nl),vsq(nl),rp(3,200)
      thkjl=thk(jl)
c
c   d1 and d2 are the x and y directions of the raypath
c
      d1=(xr-xe)/delta
      d2=(yr-ye)/delta
      nrpd=kk-jl+1
c
c   nrpd= number of raypoints in down-going part of ray path
c
      nrp=nrpd+kk
c
c   down-going part of path
c
      thk(jl)=thkjl-tkj
      rp(1,1)=xe
      rp(2,1)=ye
      rp(3,1)=ze
c   m is the number of the layer above raypoint i
c
      do 23204 i=2,nrpd
      m=jl+i-2
c
c   tng is the tangent of the incident angle in layer m
c
      tng=v(m)/sqrt(vsq(kk)-vsq(m))
      rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
      rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
      rp(3,i)=rp(3,i-1)+thk(m)
23204 continue
23205 continue
c
c   up-going part of path
c
      rp(1,nrpd+1)=rp(1,nrpd)+(delta-didjkk)*d1
      rp(2,nrpd+1)=rp(2,nrpd)+(delta-didjkk)*d2
      rp(3,nrpd+1)=rp(3,nrpd)
      thk(jl)=thkjl
      nrp1=nrp-1
c
c   nrpd2 is the number of the raypoint at the
c   top of layer (kk-1)
c
      nrpd2=nrpd+2
cuk      do 23206i=nrpd2,nrp1
      do 23206 i=nrpd2,nrp1+1  ! calculate up to the top, then compute error !!
c
c   m is the number of the layer below raypoint i
c
      m=(kk-1)-(i-nrpd2)
c
c   tng is the tangent of the incidence angle in layer m
c
      tng=v(m)/sqrt(vsq(kk)-vsq(m))
      rp(1,i)=rp(1,i-1)+thk(m)*tng*d1
      rp(2,i)=rp(2,i-1)+thk(m)*tng*d2
      rp(3,i)=rp(3,i-1)-thk(m)
23206 continue
c
      refrerr=sqrt( (rp(1,nrp)-xr)**2 + (rp(2,nrp)-yr)**2 +
     &              (rp(3,nrp)-zr)**2 ) *1000.
      rp(1,nrp)=xr
      rp(2,nrp)=yr
      rp(3,nrp)=zr
      return
      end ! of subr. refpath
c
cek    end of vel_ray.f
c
cek    begin of vel_tools.f
c
      subroutine CASEFOLD(cn) ! Urs Kradolfer, 20. Feb. 1990
c
c Input : character-string of any length
c Output: same character-string, but all letters are CAPITAL now
c         (other characters are not changed)
c
c   --> this subr. needs the function TRIMLEN !
c
c Call  : call CASEFOLD(charstring)
c
      implicit none
      character cn*(*)
      integer i, trimlen,ilen, ival
      ilen=trimlen(cn)
      do i=1,ilen
         ival=ICHAR(cn(i:i))
         if(ival.ge.97.and.ival.le.122)then
            cn(i:i)=CHAR(ival-32)
         endif
      enddo
      RETURN
      end ! of subr. casefold
c
c
      subroutine DELAZ(alat,alon,blat,blon,del,dist,az)
c
c      computes distance and azimuth from a to b
c      a and b are in decimal degrees and n-e coordinates
c      del -- delta in degrees
c      dist -- distance in km
c      az -- azimuth from a to b clockwise from north in degrees
c
      implicit none
      real alat,alon,blat,blon,del,dist,az
c
      real aa,bb,bc,olat,olon,sint,cost,rotate
      common/corect/ aa,bb,bc,olat,olon
      common/coords/ sint,cost,rotate
      COMMON/ORIGI/ REARTH,ELLIP,RLATC,RAD          !  neu
      DOUBLE PRECISION REARTH,ELLIP,RLATC,RAD       !  neu
c
cc      double precision pi2,rad,flat
      double precision pi2,   flat
      double precision alatr,alonr,blatr,blonr
      double precision tana,geoa,acol,tanb,geob,bcol
      double precision diflon,cosdel,delr,top,den,azr,colat,radius
      double precision dtan,datan,dsin,dcos,dacos,datan2
      data pi2/1.570796d0/
cc      data rad/1.745329d-02/   ! defined in subr. SETORG
calt      data flat/.993231d0/
ctest      data flat/.99330647d0/ ! neu wie in setorg
ccc      call setorg(orlat,orlon)  ! is called separately !!
      flat=rlatc !!! from subr. SETORG
c-----convert to radians
      alatr=alat*rad
      alonr=alon*rad
      blatr=blat*rad
      blonr=blon*rad
c-----convert latitudes to geocentric colatitudes
      tana=flat*dtan(alatr)
      geoa=datan(tana)
      acol=pi2-geoa
      tanb=flat*dtan(blatr)
      geob=datan(tanb)
      bcol=pi2-geob
c-----calculate delta
      diflon=blonr-alonr
      cosdel=dsin(acol)*dsin(bcol)*dcos(diflon)+dcos(acol)*
     &dcos(bcol)
      delr=dacos(cosdel)
c-----calculate azimuth from a to b
      top=dsin(diflon)
      den=(dsin(acol)/dtan(bcol))-dcos(diflon)*dcos(acol)
      azr=datan2(top,den)
c----- convert to degrees
      del=delr/rad
      az=azr/rad
      if(az.lt.0.0) az=360.+az
c-----compute distance in kilometers
      colat=pi2-(alatr+blatr)/2.d0
cold      radius=6371.277d0*
cold     & (1.d0+(3.37853d-3)*((1.d0/3.d0)-((dcos(colat))**2)))
      radius=6378.163d0*      !  neu wie in setorg
     & (1.d0+(3.35278d-3)*((1.d0/3.d0)-((dcos(colat))**2))) ! neu wie in setorg
      dist=delr*radius
      return
c
c  1/298.26 = 0.0033527795   =  flattening   ( WGS72 )
c
      end ! of subr. delaz
c
      subroutine GEOKO(x,y,xlat,xlon,i)    ! Urs Kradolfer, 24.3.87
c                         ->         1
c                         <-        -1
c
c Conversion of Swiss-coordinates X,Y [km] to LAT,LON [deg]    i = +1
c Conversion of LAT,LON [deg] to Swiss-coordinates X,Y [km]    i = -1
c
      implicit none
      real x,y,xlat,xlon
      integer i
      real seichmy
      if(.not.(i.eq.1.or.i.eq.-1)) stop'GEOKO>>> specify conversion !!!'
c
      if(i.eq.1)  call EBELL(x,y,xlon,xlat,seichMY)
      if(i.eq.-1) call ELLEB(xlon,xlat,x,y)
c
      end
c
      subroutine ELLEB(L,B,x,y)
      implicit none
      real b,x,y
      REAL L
C
C  SCHWEIZ. PROJEKTIONSSYSTEM  FORMELN VON H. ODERMATT
C  TRANSFORMATION ELLIPSOID - EBENE
C  L,B  UNTERSCHIEDE IN LAENGE UND BREITE ZU BERN (SEXAG.SEK.)
C  Y,X LANDESKOORDINATEN IN METER ( BERN 0/0)
C  MY  MERIDIANKONVERGENZ ( SEXAG. SEK.)
C
      integer i
      DOUBLE PRECISION D(5),E(5),F(5),RW(5),IW(5),P,Q,A,C,bb,bl
      DATA             BB,BL/169028.66,26782.5/
      A= L*3600.-BL
      C= B*3600.-BB
      D(1)= 0.68382546262761
      E(1)= 2.363591607471D-2
      D(2)=-3.91798328045D-8
      E(2)= 0.
      D(3)= 1.4965410352D-15
      E(3)= 4.527219881D-17
      D(4)=-8.039471422D-23
      E(4)=-3.89081120D-24
      D(5)= 7.0021390D-30
      E(5)= 2.3407700D-31
      F(1)= 4.515344386039D1
      F(2)= 1.17912305209D-4
      F(3)= 5.8474201864D-10
      F(4)= 2.73386187D-15
      F(5)= 1.4308547D-20
C
      P=30.91849390613 * A
      Q=C*F(5)
      I=5
  1   I=I-1
      Q=C*(Q+F(I))
      IF(I.GT.1) GOTO 1
C
      RW(1)=Q
      IW(1)=P
      DO 2 I=2,5
      RW(I)=Q*RW(I-1)-P*IW(I-1)
      IW(I)=P*RW(I-1)+Q*IW(I-1)
  2   CONTINUE
C
      y=D(5)*RW(5)
      x=D(5)*IW(5)
      I=5
  3   I=I-1
      y=y+D(I)*RW(I)
      x=x+D(I)*IW(I)
      IF(I.GT.1) GOTO 3
      y= y/1000.+200.       !  output in kilometers
      x= x/1000.+600.
C
      RETURN
      END ! of subr. ELLEB
c
      subroutine EBELL(YL,XB,L,B,MY)
C
C  SCHWEIZ. PROJEKTIONSSYSTEM FORMELN VON H. ODERMATT
C  TRANSFORMATION EBENE ELLIPSOID
C  YL,XB LANDESKOORDINATEN IN METER (BERN 600/200)
C  L,B LAENGE UND BREITE  (GRAD)
C  MY MERIDIEANKONVERGENZ (SEXAG.SEK.)
C
      implicit none
      REAL YL,XB,L,B,MY
      integer i
      DOUBLE PRECISION A(8),BB(8),C(8),RZ(8),QZ(8),P,Q,X,Y,DB,DL
cccc      dimension A(8),BB(8),C(8),RZ(8),QZ(8)
cccc      real P,Q,X,Y,DB,DL
C
      X=1000.*XB-200000.
      Y=1000.*YL-600000.
C
      A(1)=1.4623614572021
      BB(1)=3.4564252673326D-2
      A(2)=1.225255821052D-7
      BB(2)=2.89600437564D-9
      A(3)=1.3687923002D-14
      BB(3)=4.651046030D-16
      A(4)=1.971224191D-21
      BB(4)=6.43850954D-23
      A(5)=2.97898051D-28
      BB(5)=9.600412D-30
      A(6)=4.650273D-35
      BB(6)=1.50512D-36
C     A(7)=7.48203D-42
      A(7)=0.D00
C     BB(7)=2.422D-43
      BB(7)=0.D00
C     A(8)=1.229D-48
      A(8)=0.D00
      C(1)= 2.2146704979846D-2
      C(2)=-1.280815253730D-9
      C(3)= 7.4775676024D-18
      C(4)= 4.691943327D-24
      C(5)=-3.6550101D-31
      C(6)= 3.71615D-39
C     C(7)= 1.6901D-45
      C(7)= 0.D00
C     C(8)= 1.96D-52
      C(8)= 0.D00
      RZ(1)= X
      QZ(1)= Y
cuk      DO 1 I=2,8   !  sufficiant to loop until 6 only
      do 1 I=2,6      !  otherwise VAX has an overflow...
      RZ(I)= X*RZ(I-1)-Y*QZ(I-1)
      QZ(I)= Y*RZ(I-1)+X*QZ(I-1)
    1 CONTINUE
      I=8
      Q=A(I)*RZ(I)
      P=A(I)*QZ(I)
      MY=0.
    2 I=I-1
      Q=Q+A(I)*RZ(I)
      P=P+A(I)*QZ(I)
      MY=MY+BB(I)*QZ(I)
      IF(I.GT.1) GOTO 2
      DL=3.2343101932327D-2 * P
      DB=Q*C(8)
      I=8
    3 I=I-1
      DB=Q*(DB+C(I))
      IF(I.GT.1) GOTO 3
C
cc      write(6,*)'vor L= und B= in subr. ebell!'
      L=(DL + 26782.5D00)/3600.
      B=(DB + 169028.66D00)/3600.
      RETURN
      END ! of subr. EBELL
c
      subroutine SDC(x,y,xlat,xlon,i)     ! 7. Jan. 1988
c
c     SDC needs subroutines SETORG, DIST and REDIST
c
c
      implicit none
      real x,y,xlat,xlon
      integer i
c
      real olat,olon,aa,bb,bc,sint,cost,rotate
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
c
cc      call setorg()  ! is called separately !!
c
      if(.not.(i.eq.-1.or.i.eq.1)) stop'SDC>>> specify conversion !!!'
c
      if(i.eq.1)  call redist(x,y,xlat,xlon)
      if(i.eq.-1) call dist(xlat,xlon,x,y)
c
      end ! of subr. SDC
c
      subroutine DIST(xlat,xlon,xkm,ykm)   ! 7. Jan. 1988
c
      implicit none
      real xlat,xlon,xkm,ykm
c      
      real olat,olon,aa,bb,bc,sint,cost,rotate
      real q,yp,xx
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
      DOUBLE PRECISION LAT1,LAT2,LAT3
c     conversion of latitude and longitude to kilometers relative
c     to center of coordinates
c
      q=60*xlat-olat
      yp=q+olat
      lat1=datan(rlatc*dtan(RAD*yp/60.0))
      lat2=datan(rlatc*dtan(RAD*OLAT/60.0))
      LAT3=(LAT2+LAT1)/2.
      xx=60*xlon-olon
      q=q*aa
      xx=xx*bb*dcos(LAT3)
cc** rotate coordinate system clockwise
      yp=cost*q-sint*xx
      xx=cost*xx+sint*q
      q=yp
c
      xkm=xx
      ykm=q
c
      return
      end ! of subr. dist
c
      subroutine REDIST(Xkm,Ykm,XLAT,XLON)   ! 7. Jan. 1988
C   by Edi
C CALCULATES DISTANCE BETWEEN TWO POINTS,  INPUT IN RECTANGULAR XY-COOR.
C   OUTPUT IN LATITUDE/LONGITUDE
C
C     CONVERSION OF  KILOMETERS TO LATITUDE AND LONGITUDE
c
      implicit none
      real xkm,ykm,xlat,xlon
c
      real olat,olon,aa,bb,bc,sint,cost,rotate
      real xx,yy,y,x,q,lat,yp,bcl,p,lon
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
      DOUBLE PRECISION LAT1,LAT2,LAT3,clat1       ! neu
C
c
      xx=xkm
      yy=ykm
c*****	Following code is commented-out
c	and replaced by correct code to undo rotation
c	WLE 1/3/95

c      rphi=0.0   ! 5.4.91 u.k.
cc
c      SINT=DSIN(RPHI*RAD)
c      COST=DCOS(RPHI*RAD)
c      TANT=DTAN(RPHI*RAD)
cC---- ROTATE COORDINATES CLOCKWISE BACK BY RPHI
c      IF(RPHI.LT.0.1) GOTO 120
c      Y=YY*COST-XX*SINT
c      X=XX/SINT+Y*TANT
c      GOTO 130
c  120 X=XX
c      Y=YY

c	rotate coordinate system into original orientation

	x=xx*cost-yy*sint
	y=yy*cost+xx*sint
c	WLE 1/3/95

  130 CONTINUE
      IF(ABS(AA).LT.0.0000001) GOTO 900
      Q=Y/AA
      LAT=(Q+OLAT)/60.
      XLAT=Q+OLAT - 60.*LAT
      YP=60.*LAT+XLAT
      LAT1=DATAN(RLATC*DTAN(YP*RAD/60.0))
      LAT2=DATAN(RLATC*DTAN(OLAT*RAD/60.))
      LAT3=(LAT1+LAT2)/2.
      CLAT1=DCOS(LAT3)
      BCL=BB*CLAT1
      IF(ABS(BCL).LT.0.000001) GOTO 900
      P=X/(BB*CLAT1)
      LON=(P+OLON)/60.
      XLON=P+OLON - 60.*LON
      xlat=lat+xlat/60.
      xlon=lon+xlon/60.
      RETURN
  900 WRITE(6,1000) AA,BB,CLAT1
 1000 FORMAT(/,2X,' SUBR. REDISt: AA=',F10.5,2X,'BB=',F10.5,2X,
     1'COS(LAT1)=',F10.7,5X,'DIVISION BY ZERO, RUN STOPS HERE',/)
      stop'REDIST>>> division by zero!!'
      END ! of subr. redist
c
      subroutine SETORG(orlat,orlon,rrotate,ifil) 
c
      implicit none
      real orlat,orlon,rrotate
      integer ifil
c
      real olat,olon,aa,bb,bc,sint,cost,rotate
      integer icoordsystem
      common/GEO_COORSYSTEM/olat,olon,
     &                      rearth,ellip,rlatc,rad,
     &                      aa,bb,bc,sint,cost,rotate, icoordsystem
      double precision rearth,ellip,rlatc,rad
c
c
      double precision r,lat1,lat2,dela,delb,PHI,BETA
c
c     O(r)LAT & O(r)LON : origin of cartesian coordinate system
c      north  w e s t
c
c
cEK Dez94      rotate=0.0  !......
c
      rotate=rrotate
cEK
      if(orlat.eq.0.0.and.orlon.eq.0.0)then
         olat=46.95240  ! BERN North
         olon=-7.439583  ! BERN West
      else
         olat=orlat
         olon=orlon
      endif
c
      olat=olat*60. ! minutes N
      olon=olon*60. ! minutes W
c
       rad=0.017453292d0
C
C  OLD ELLIPSOID:
C      rlatc = 0.99330647d0
C      data rearth /6378.163d0/, ellip /298.26d0/
C
C  NEW ELLIPSOID FOR WHOLE EARTH:   WGS72 == WORLD GEODETIC SYSTEM 1972
C
C  ALSO SET RLATC ACCORDING TO ORIGIN
C
      REARTH=6378.135D0
      ELLIP=298.26         ! (flattening)
C
c rlatc = tan(geocentric LAT) / TAN(GEOGRAPHIC LAT)
C
C
C CALCULATE RLATC:  CONVERSION FROM GEOGRAPHICAL LAT TO GEOCENTRICAL LAT
C
      PHI=OLAT*RAD/60.0              !  phi=geogr. lat
      BETA=PHI-DSIN(PHI*2.)/ELLIP    !  beta=geoc. lat
      RLATC=DTAN(BETA)/DTAN(PHI)
c
C
C WRITE ELLIPSOIDE CONSTANTS
C
      if(ifil.gt.0)then
         write(ifil,*)
         write(ifil,*)
         write(ifil,*)'SHORT DISTANCE CONVERSION on ELLIPSOIDE of'//
     &                ' WORLD GEODETIC SYSTEM 1972 (WGS72)'
         write(ifil,*)'=========================================='//
     &                '==================================='
         write(ifil,*)
         write(ifil,'('' Radius at equator (REARTH)= '',f10.5,
     &                ''  km'')') rearth
         write(ifil,'(''   1. / (ellipticity)      = '',f10.3)') ellip
         write(ifil,*)
         write(ifil,*)'Origin of cartesian coordinates [degrees]:'
         if(orlat.eq.0.0.and.orlon.eq.0.0)then
            write(ifil,*)' (Origin = city of BERNE, Switzerland)'
         endif
         write(ifil,'(1x,f12.7,'' N'',5x,f12.7,'' W'')')
     &               olat/60.,olon/60.
         write(ifil,*)
         write(ifil,*)' Rotation angle (in degr.) clockwise from'
c        write(ifil,'(''North  rotate= '',f6.1,)') rotate
         write(ifil,*) rotate
         write(ifil,*)
      endif
c
c   calculate aa &  bb
c   length of one minute of lat and lon in km
c
      lat1=datan(rlatc*dtan(olat*rad/60.0))       ! geoc. lat for OLAT
      lat2=datan(rlatc*dtan((olat+1.)*rad/60.0))  ! geoc. lat for (OLAT+1min)
      dela=lat2 - lat1
      r=rearth*(1.0 - (dsin(lat1)**2)/ellip)      ! kugelradius fuer lat=OLAT
      aa=r*dela   ! aa = 1 min geogr. lat
      delb=dacos(dsin(lat1)**2 + dcos(rad/60.)*dcos(lat1)**2)
      bc=r*delb     ! bc = 1 min geogr. lon
      bb=r*delb/dcos(lat1)
      if(ifil.gt.0)then
         write(ifil,'('' Radius of sphere at OLAT = '',f10.3,'' km'')')r
         write(ifil,*)
         write(ifil,*)'Conversion of GEOGRAPHICAL LATITUDE to '//
     &                'GEOCENTRICAL LATITUDE:'
         write(ifil,*)'RLATC = TAN(GEOCENTR.LAT) / TAN(GEOGRAPH.LAT)'
         write(ifil,'(1x,''RLATC = '',f12.8)') rlatc
         write(ifil,*)
         write(ifil,4) aa, bc
 4       format (10x,'Short distance conversions:',/,
     &           10x,'one min lat ~ ', f7.4,' km ',/,
     &           10x,'one min lon ~ ', f7.4,' km ',/)
         write(ifil,*)
         write(ifil,*)
      endif
C
c***  convert coordinates with rotation cosines
      sint=sin(rotate*rad)
      cost=cos(rotate*rad)
31    continue
c
      return
      end ! of subr. setorg
c
      subroutine TIMECLEAR(iyr,imo,iday,ihr,imin,sec,itime)
c
c     Urs Kradolfer, Winter 1987/88
c
      implicit none
c
      integer iyr,iyr1,imo,iday,ihr,imin, juliam,itime
      real sec,sec1
c
      if(iyr.lt.1900)then
         iyr1=iyr+1900
      else
         iyr1=iyr
      endif
      itime=JULIAM(iyr1,imo,iday,ihr,imin)
   1  sec1=sec
      if(sec.lt.0)then
         sec=sec1+60.
         itime=itime-1
      endif
      if(sec.gt.60.)then
         sec=sec1-60.
         itime=itime+1
      endif
      if(sec.lt.0.0.or.sec.gt.60.0) goto 1
      call DATUM(itime,iyr1,imo,iday,ihr,imin)
      if(iyr.lt.1900)then
         iyr=iyr1-1900    ! if input-year was of form  87  give it out so !
      else
         iyr=iyr1
      endif
c
      return
c
      end ! of subr. timeclear
c
      INTEGER FUNCTION JULIAM(IYR,IMO,IDY,IHR,IMN)
C UMRECHNEN VON JAHR-MONAT-TAG-STUNDEN-MINUTEN IN MINUTEN:
C   (WENN IMN 4-BYTE INTEGER, DANN JAHR < 4000)
      implicit none
      integer iyr,imo,idy,ihr,imn
      integer KMO(12)
      integer leap,ky,km,kd,ky4,ky1,ky0,kl,l
      DATA KMO/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA LEAP/1/
      KY= IYR
      KM= IMO
      KD= IDY
      IF(KM.LE.0) KM= 1
10    JULIAM= 365*KY
      KD= KMO(KM)+KD
      KY4= KY/4
      KY1= KY/100
      KY0= KY/1000
      KL= LEAP*(KY4-KY1+KY0)
      L=0
      IF(KY4*4.EQ.KY.AND.(KY1*100.NE.KY.OR.KY0*1000.EQ.KY))L= LEAP
      IF(L.NE.0.AND.KM.LT.3) KL= KL-LEAP
      JULIAM= JULIAM+KD+KL
      JULIAM= JULIAM*24+IHR
      JULIAM= JULIAM*60+IMN
      return
      END ! of integer function juliam
c
      subroutine DATUM(ITF,IYR,IMO,IDY,IHR,IMN)
C UMRECHNEN DES DATUMS IN MINUTEN (CF. JULIAM) IN YR-MO-DY-HR-MI
C   (MIT IMN<2**31, JAHR < 4000
c
      implicit none
      integer iyr,imo,idy,ihr,imn
c      
      integer id,l,iyr4,iyrh,iyrt,ld,i
      integer KMO(12)
      INTEGER ITF,K,KH
      DATA KMO/31,28,31,30,31,30,31,31,30,31,30,31/
      K= ITF/60
      IMN= ITF-K*60
      KH= K/24
      IHR= K-KH*24
      IYR= KH/365
5     ID= KH-IYR*365
      L= 0
      IYR4= IYR/4
      IYRH= IYR/100
      IYRT= IYR/1000
      LD= IYR4-IYRH+IYRT
      IF(IYR4*4.EQ.IYR.AND.(IYRH*100.NE.IYR.OR.IYRT*1000.EQ.IYR)) L= 1
      ID= ID-LD+L
      IF(ID.GT.0) GOTO 10
      if(id.eq.0.and.ihr.eq.0.and.imn.eq.0) then
          idy= 0
          imo= 0
          return
      endif
      IYR= IYR-1
      GOTO 5
10    KMO(2)= 28+L
      DO 20 I=1,12
      ID= ID- KMO(I)
      IF(ID.LE.0) GOTO 30
20    CONTINUE
      I=12
30    IDY= ID+KMO(I)
      IMO= I
      RETURN
      end ! of subr. datum
c
      integer function TRIMLEN(t)   ! Urs Kradolfer, June 1986
c
c     Call:    nc=TRIMLEN(char)
c
c          --> nc says, how many characters the input-string has
c              (ignoring trailing blanks!).
c
      implicit none
      character t*(*)
      do 1 trimlen=LEN(t),1,-1
    1    if(t(trimlen:trimlen).ne.' ')RETURN
      trimlen=1
      end ! of integer function trimlen
c
      subroutine FREEUNIT(iunit)
      implicit none
c
c     Purpose: Get a free FORTRAN-unit
c
c     Output: iunit  a free (file-)unit number
c
c     Urs Kradolfer, 16.7.90
c
      integer iunit
      logical lopen
c
      do iunit=10,999
         if(iunit.eq.999) stop'FREEUNIT>>> no free unit found!'
         inquire(unit=iunit,opened=lopen)
         if(.not.lopen) RETURN
      enddo
c
      RETURN
c
      end ! of subr. freeunit
c
cek    end of vel_tools.f
c
cek    begin of vel_topo.f
c
      subroutine CHTOP(xx,yy,zk,
     &                 topo1file,topo2file)   !  Dummy routine
cEK
c     this routine for a given pair of coordinates (xx,yy) provides
c     the elevation (zk in km) by reading it off from topo-arrays
cEK
c     If necessary, replace it with appropriate routine and 
c     topo information
c
      implicit none
      real xx,yy,zk
      character*(*) topo1file,topo2file
c
c
      zk=0.0
      r e t u r n
      end
c
      subroutine REGION (ityp,x,y,cname,nreg,
     &                   regnamfile,regkoordfile) ! Urs Kradolfer, 6.7.87
c
c     originally from Manfred Baer
c
c     modified 6.7.87, 28.6.91 Urs Kradolfer 
c
      implicit none
c
      character*(*) regnamfile, regkoordfile
      character*32 cname, place
      integer ityp,nreg, iii, iregread
      real x,y, xlat,xlon
c
      save iregread
c
      if(iregread.eq.0)then
         call REGREAD(regnamfile,regkoordfile)
         iregread=1
      endif
c
      nreg=0
c
      if(ityp.eq.1) goto 1
      if(ityp.eq.2.or.ityp.eq.3) goto 2
      stop'REGION>>> illegal coordinate-type!'
    1 continue
      place=' '
      call REGCH(x,y,place,nreg)
      if(nreg.eq.0)then
         call GEOKO(x,y,xlat,xlon,1)
         goto 22
      else
         cname=place
      endif
      goto 999
c
    2 continue
      xlat=x
      xlon=y
  22  place=' '
c     call REGNU(0,xlat,xlon,GELAT,GELON,PLACE,NREG)
      if(ityp.eq.1) iii=0
      if(ityp.eq.2) iii=0
      if(ityp.eq.3) iii=1
      call REGWORLD(iii,xlat,xlon,place,nreg)
      if(place.eq.' ')then
         cname='*****'
      else
         cname=place
      endif
      goto 999
c
  999 continue
      return
      end ! of subr. region
c
      subroutine REGREAD(regnamfile,regkoordfile)
      implicit none
      integer indfe,lt50,lt25
      real xnrlon
      character*(*) regnamfile,regkoordfile
      character irname*24300
      common /regioncom/ indfe(730),lt50(111),lt25(421),
     &                   irname,xnrlon(14400)
c
      integer trimlen,nstart,nr
      integer nc
      integer i,j,k
      integer iu1, iu2
      character*40 cregion
cu      character*64 inputfile
c
cu      inputfile='/users/kradi/velest/regionsnamen.dat'
      call FREEUNIT(iu1)
cVMS      open(iu1,file=regnamfile,status='OLD',readonly)
      open(iu1,file=regnamfile,status='OLD')
      nstart=0
      do i=1,729
          indfe(i)=nstart
          read(iu1,'(i5,2x,a)') nr,cregion
          if(nr.ne.i) then
              write(6,'(''sequence error at:'',i5,a,i5)') nr,cregion,i
              stop'REGREAD>>> error reading file REGIONSNAMEN.DAT'
          endif
          nc=TRIMLEN(cregion)
          irname(nstart+1:nstart+nc)=cregion(1:nc)
          nstart=nstart+nc
      enddo
      indfe(730)=nstart
      do i=1,420
          lt25(i)=nstart
          read(iu1,'(i5,2x,a)') nr,cregion
          if(i.le.400) then
             k=i+999
          else
             k=1019+(i-401)*20
          endif
          if(k.ne.nr) then
              write(6,'(''sequence error at:'',i5,a,i5)') nr,cregion,i
              stop'REGREAD>>> error reading file REGIONSNAMEN.DAT'
          endif
          nc=TRIMLEN(cregion)
          irname(nstart+1:nstart+nc)=cregion(1:nc)
          nstart=nstart+nc
      enddo
      lt25(421)=nstart
      do i=1,110
          lt50(i)=nstart
          read(iu1,'(i5,2x,a)') nr,cregion
          if(i.le.100) then
              k=i+199
          else
              k=209+(i-101)*10
          endif
          if(k.ne.nr) then
              write(6,'(''sequence error at:'',i5,a,i5)') nr,cregion,i
              stop'REGREAD>>> error reading file REGIONSNAMEN.DAT'
          endif
          nc=TRIMLEN(cregion)
          irname(nstart+1:nstart+nc)=cregion(1:nc)
          nstart=nstart+nc
      enddo
      lt50(111)=nstart
      close(iu1)
c      write(6,'(''anzahl worte fuer regions namen:'',i6)') nstart
c  jetzt: flinn-engdahl-regionen
cu      inputfile='/users/kradi/velest/regionskoord.dat'
      call FREEUNIT(iu2)
cVMS      open(iu2,file=regkoordfile,status='OLD',readonly)
      open(iu2,file=regkoordfile,status='OLD')
      do j=1,14400,8
          read(iu2,'(8f9.3)') (xnrlon(i),i=j,j+7)
      enddo
      close(iu2)
c
      return
      end ! of subr. regread
c
      subroutine REGCH(x,y,place,nreg)
      implicit none
      real x,y
      integer nreg
      integer indfe,lt50,lt25,nx,ny,nrpu,i1,i2
      real xnrlon,dx,dy,xx,yy
      character irname*24300
      COMMON /REGIONcom/ INDFE(730),LT50(111),LT25(421),
     &                   irname,XNRLON(14400)
      character place*32
c
      nreg=0
      IF(X.LT.480.) RETURN
      IF(X.GE.865.) RETURN
      IF(Y.LE. 62.) RETURN
      IF(Y.GT.302.) RETURN
      dx=17.5
      dy=12.0
      xx=X-480.
      yy=302.-Y
      nx=xx/dx
      NY=yy/dy
      nrpu=ny*20+nx
      nreg=nrpu
      if(nx.ge.20) then
        nrpu=400+ny
        nreg=ny*20 + 19
      endif
      i1=lt25(nrpu+1)+1
      i2=lt25(nrpu+2)
      nreg=nreg+1000
      place=irname(i1:i2)
      return
      END ! of subr. regch
c
      subroutine REGWORLD(ityp,cord1,cord2,place,nreg)
      implicit none
      integer ityp,nreg,last,num,in,nlo,latin,mhemis,i,iuk
      integer nval,i1,i2
      real cord1,cord2,xnrlon,galat,galon,pi
      real galor,galar,gelat,gelon,xazr,xdr,q,as
      real pirad,pi2,gelac,xlatr,gelor,ck1,ck2,ac
      real alon,xlat,xlon,xreg
      integer indfe,lt50,lt25
      character irname*24300
      COMMON /REGIONcom/ INDFE(730),LT50(111),LT25(421),
     &                   irname,XNRLON(14400)
C     PROGRAM TO COMPUTE REGION NUMBER AND NAME BASED ON THE INPUT
C        ITYP <=0         LATITUDE&LONGITUDE
C        ITYP > 0         DISTANCE & AZIMUTH
C    THE REGION NAMES ARE FROM FLINN & ENGDAHL#S TABLES
      character place*32
      real BUF(40)
      character*1 cbk
      DATA GALAT/46.7706/,GALON/8.09428/
      data cbk/' '/
C
      PI=3.141593
      PIRAD=PI/180.
      PI2=PI*2.
      LAST=89
C
C     CONVERT ARRAY COORDINATES TO RADIANS
      GALOR=GALON*PIRAD
      GALAR=GALAT*PIRAD
      IF(ITYP) 850,850,860
C
  850 GELAT=CORD1
      GELON=CORD2
      NUM=1
      GO TO 212
C
  860 XAZR=CORD2*PIRAD
      XDR= CORD1*PIRAD
      NUM=1
C
      Q=SIN(GALAR)*COS(XDR)+COS(GALAR)*SIN(XDR)*COS(XAZR)
      IF(ABS(Q).GT.1.) GO TO 891
      AS=ASIN(Q)
      GO TO 890
C
  891 CONTINUE
      GELON=0.
      GELAT=0.
      NUM=NUM-1
      GO TO 212
C
  890 GELAC=AS
      XLATR=GELAC
      IF(COS(GELAC))45,46,45
   45 IF(COS(GALAR)) 47,46,47
   46 GELOR=0.0
      GO TO 34
C
   47 CONTINUE
      CK1=(COS(XDR)-SIN(GALAR)*SIN(XLATR))/
     *(COS(GALAR)*COS(XLATR))
C
      CK2=(SIN(XAZR)*SIN(XDR))/COS(XLATR)
      IF(CK1.GE.0.) GO TO 37
      IF(CK2.GE.0.) GO TO 882
C
      IN=1
      GO TO 881
C
  882 IN=2
      GO TO 883
C
   37 IF(CK2.GE.0.) GO TO 885
C
      IN=3
      GO TO 881
C
  885 IN=4
C
  883 IF(ABS(CK1).GT.1.0) GO TO 887
      AC=ACOS(CK1)
      GO TO 888
C
  881 IF(ABS(CK2).GT.1.0) GO TO 887
      AS=ASIN(CK2)
C
  888 GO TO (33,32,31,36),IN
C
   33 GELOR=GALOR-PI-AS
      GO TO 34
C
   32 GELOR=GALOR+AC
      GO TO 34
C
   31 GELOR=GALOR+AS
      GO TO 34
C
   36 GELOR=GALOR+AC
      GO TO 34
C
  887 CONTINUE
      GELOR=0.0
C
C     CONVERT LONGITUDE AND LATITUDE FROM RADIANS TO DEGREES
C     AND REDUCE LONGITUDE TO A FORM LESS THAN 360 DEGREES
C
   34 CONTINUE
      ALON=ABS(GELOR)
  506 IF(ALON.LT.PI2) GO TO 116
      ALON=ALON-PI2
      GO TO 506
C
  116 GELOR=SIGN(ALON,GELOR)
      IF(ABS(GELOR).GT.PI) GELOR=GELOR-SIGN(PI2,GELOR)
C
C  LONGITUDE IN DEGREES
C
      GELON=GELOR/PIRAD
C
C     LATITUDE IS CONVERTED FROM GEOCENTRIC TO GEODETIC COORDINATES
C
      GELAT=GELAC/PIRAD
C
  212 IF(NUM.EQ.1) GO TO 5000
      place=' '
      GO TO 900
C
 5000 CONTINUE
C
C      CLEAR BUFFER
C
C     ASSIGN REGION NUMBER AND NAME
C
      place=' '
  215 XLAT=GELAT
      XLON=GELON
      IF(XLON.LT.0.) XLON=XLON+360.
C
C     DETERMINE EARTH QUADRANT
C
      IF(XLAT.LT.0.) GO TO 4
      IF(XLON.GT.180.) GO TO 6
C
C     LATITUDE POSITIVE,LONGITUDE 0-180
C
      MHEMIS=0
      GO TO 9
C
C     LATITUDE POSITIVE,LONGITUDE 181-360
C
    6 MHEMIS=3600
      XLON=360.-XLON
      GO TO 9
C
    4 IF(XLON.GT.180.) GO TO 8
C
C     LATITUDE NEGATIVE,LONGITUDE 0-180
C
      MHEMIS=7200
      GO TO 9
C
C     LATITUDE NEGATIVE,LONGITUDE 181-360
C
    8 MHEMIS=10800
      XLON=360.-XLON
C
    9 CONTINUE
C
C     CHECK FOR ZERO LATITUDE
C
      IF(ABS(XLAT).GE.1.0) GO TO 13
      LATIN=0
      GO TO 17
C
   13 LATIN=ABS(XLAT)
      IF(LATIN.GT.LAST) LATIN=LAST
C
   17 CONTINUE
C
      NLO=ABS(XLON)
      I=MHEMIS + 40*LATIN + 1
      do iuk=1,40
         buf(iuk)=xnrlon(i+iuk-1)
      enddo
      NVAL=BUF(1)*1000.+1.1
      DO 56 I=2,NVAL
      IF(NLO.LT.int(BUF(I))) GOTO 57
56    CONTINUE
      I=NVAL+1
57    NREG=BUF(I-1)
      XREG=BUF(I-1)-NREG
      NREG=XREG*1000. + 0.1
      I1=INDFE(NREG)+1
      I2=INDFE(NREG+1)
      place=irname(i1:i2)
C
  900 continue
      return
C
      END ! of subr. regworld
c
cek    end of vel_topo.f
c
cek    begin of matrinv.f
c
      SUBROUTINE matrinv(N,A,B)
c
c MATRINV invertiert eine SYMMETRISCHE NxN-Matrix
c
c Aufruf: call MATRINV(N,A,B)
c
c         N = Anz. Zeilen/Kolonnen der Matrix
c         A = zu invertierende Matrix (Input)
c         B = inv(A)    (Output)
c
      DIMENSION A(N*N),B(N*N)
      B(1)= 1.0/A(1)
      IF(N-1) 60,60,2
2     NN= N*N
      DO 4 I=2,NN
4     B(I)= 0.0
      MM=1
      KN= 0
      DO 50 M=2,N
      K= M-1
      MM= MM+N+1
      KN= KN+N
      EK= A(MM)
      MI= M-N
      DO 10 I=1,K
      MI= MI+N
      IJ= I-N
      JM= KN
      DO 10 J= 1,K
      IJ= IJ+N
      JM= JM+1
10    EK= EK-A(MI)*B(IJ)*A(JM)
      B(MM)= 1.0/EK
      MI= M-N
      IM= KN
      DO 30 I=1,K
      IM= IM+1
      IJ= I-N
      JM= KN
      DO 20 J= 1,K
      IJ= IJ+N
      JM= JM+1
20    B(IM)= B(IM)-B(IJ)*A(JM)*B(MM)
      MI= MI+N
30    B(MI)= B(IM)
      IM= KN
      DO 40 I=1,K
      IM= IM+1
      MJ= M-N
      IJ= I-N
      DO 40 J=1,K
      MJ=MJ+N
      IJ= IJ+N
40    B(IJ)= B(IJ)+B(IM)*B(MJ)*EK
50    CONTINUE
60    continue
      return
      END
c
cek    end of matrinv.f
c

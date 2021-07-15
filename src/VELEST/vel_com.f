c     file VELESTCOMMON.FOR                        Version ETH-19SEP91
c******************************************
c     common-blocks of program VELEST     *  by Urs Kradolfer
c******************************************
c
c
c-----------------------------------------------------------------------------
c
c     IN THIS SECTION YOU CAN CHANGE THE ARRAY-SIZES:
c
      integer ieq,inshot,itotmodels,inltot,ist,
     &        maxobsperevent 
      parameter(ieq=5000)      ! number of earthquakes (Simultaneous mode)
c      parameter(ieq=1)        ! number of earthquakes (single event mode ONLY)
      parameter(inshot=50,    ! max. number of shots (for simultaneous mode)
     &          itotmodels=2, ! max. number of velocity models
     &          inltot=100,    ! max. number of layers per velocity model
     &          ist=650,      ! max. number of stations in stationlist
     &          maxobsperevent=180) ! max. number of observations per event
c
c
c-----------------------------------------------------------------------------
c
      integer ilip,inva,ikva,inrpmax
      parameter(ilip=(4*ieq+inshot+inltot+ist))
      parameter(inva=(ilip-1))
      parameter(ikva=(inva*ilip/2))
      parameter(inrpmax=(2*inltot+2))
c
c     Control-parameters:
c
      integer nvar,nvareff,kvar,nsta,ksta,legs,lip,neqs,
     &        nsinv,nshcor,nshfix,icount,nshot 
      common/DIMEN/nvar,nvareff,kvar,nsta,ksta,legs,lip,neqs,nshot,
     &             nsinv,nshcor,nshfix, icount
      integer ittmax,invertratio,iresolcalc,ifixsolution
      real    zmin,zmininput,delmin,veladj,swtfac,xythet,
     &        zadj,vthet,othet,zthet,stathet,rmsmin,dmax 
      common/PARAM/ zmin,zmininput,ittmax,delmin,veladj,swtfac,xythet,
     &              invertratio,iresolcalc,
     &              zadj,vthet,othet,zthet,stathet,rmsmin,
     &              dmax,ifixsolution
      integer icnvout,istaout,ismpout,irayout,idrvout,ialeout,idspout,
     &        irflout,irfrout,iresout 
      common/OUTPUTFILES/ icnvout,istaout,ismpout,
     &                    irayout,idrvout,ialeout,idspout,
     &                    irflout,irfrout,iresout
      real zshift
      common/CORR/ zshift
      logical single_turbo
      integer icoordsystem,itrial,ised,isingle,itopo,iturbo,
     &        lowveloclay,ielev,iuseelev,iusestacorr
      real    ztrial
      common/COORDSYST/ icoordsystem, itrial, ztrial, ised, isingle,
     &                  itopo, iturbo,lowveloclay, single_turbo,
     &                  ielev(ist), iuseelev, iusestacorr
c
c     Model-parameters:
c
      integer nmod,nsp,nltot,nplay,laysum
      common/MODELA/ nmod,nsp,nltot,nplay(inltot),laysum(inltot)
      real    v,vsq,h,thk,tkj,delta
      integer jl
      common/LAYER/ v(inltot),vsq(inltot),h(inltot),thk(inltot),jl,tkj,
     &              delta
      real    vp,hp,vdamp,
     &        thkp,vpvs
      integer ireflector,lmax
      character*1 reflchar
      common/LAYVEL/ vp(itotmodels,inltot),hp(itotmodels,inltot),
     &               vdamp(itotmodels,inltot),
     &               thkp(itotmodels,inltot),vpvs,
     &               ireflector,lmax, reflchar       ! lmax added 5.4.91 uk
c
c     Station-parameters:
c
      real     ptcor,stcor,d,x
      integer  map1,map2,model 
      common/STATN/ ptcor(ist),stcor(ist),d(ist,3,ieq),
     &              model(itotmodels*ist),    ! to be verified !!! (hrm 30.3.92)
     &              x(ist,3),map1(ist),map2(ist)
      real     xla,xlo
      common/STLTLN/ xla(ist),xlo(ist)
cccc      character*1 cns,cew
      character*1 prmk
      real        pt,w,sphase,
     &            depthsofinput,tcalc,amx,prx
      integer     kpwt,istm,iphase,knobs,
     &            iain,nobswithw0,idelta,nactualsta,
     &            nstaeff
      common/OBB/ pt(ist,ieq),kpwt(ist,ieq),istm(ist,ieq),w(ist,ieq),
     &            iphase(ist,ieq),sphase(ist,ieq),knobs(ieq),
     &            depthsofinput(ieq),
     &            tcalc(ist),iain(ist),
     &            nobswithw0,amx(ist),prx(ist),idelta(ist),
     &            prmk(ist,2),nactualsta(ist),nstaeff
      real tctime
      common/EXTRA/ tctime(ist,ieq)
c89      character*4 smn,stn,blank,blank0
      character*6 smn,stn,blank,blank0
c89      data blank,blank0 /'    ','0   '/
      common/CHARc/ smn(ist,ieq),stn(ist),blank,blank0
c
c     Event-parameters:
c
      character*32 regionname
      integer      iyr,imo,iday,ihr,imin,ifx,isconstrain,iconstrain,
     &             nmag,nreg
      real         e,emag,xmagni,xmagnitude,sdxmagnitude
      common/EVENT/iyr(ieq),imo(ieq),iday(ieq),ihr(ieq),imin(ieq),
     &             e(4,ieq),emag(ieq),ifx(ieq),
     &             xmagni(ist),nmag,xmagnitude,sdxmagnitude,
     &             regionname,nreg,
     &             isconstrain(3), iconstrain(ieq)
ccc     &             , effdeltaz(ieq)
      integer igap
      common/GAP/ igap(ieq)   ! max. ieq events
c
c     'Normal-equations-parameters':
c
      real rht,res,rms,b,s,avres,steplen
      character*80 smpline
      common/A4/ rht(inva*2),res(ist,ieq),rms(ieq),b(inva*2),s(inva*2),
     &           avres(ieq), steplen,smpline
      real  scale
      common /FIX1/scale(7)
      real g
      common /GMATRIX/g(ikva*2)
      real gcopy
      common /GSAVE/gcopy(60,inva*2)
      real dtdv,dtdr
      common/DTTX/ dtdv(inltot),dtdr(3)
      real gg2,ggt,ggg,gtg,ggti,drm
      common/DATARESOL/gg2(maxobsperevent,4),ggt(4,maxobsperevent),
     &                 ggg(4,maxobsperevent),gtg(4,4),ggti(4,4),
     &                 drm(maxobsperevent,maxobsperevent)
c
      real Rc,Rs,SV,COVc,COVs,ale,R3,GG
      common/goodness/ Rc(4,4),Rs(4,4), SV(4), COVc(4,4),COVs(4),
     &                 ale(ieq),R3(3,3),
     &                 GG(100,100)
c
c     Statistics-parameters:
c
      integer irefrlayer,ihypoclayer,noheadwave,irefllayer
      real    refraylen,avhraylen,avvraylen,hitlay,sterr,
     & direrr,refrerr,reflerr
      common /HITL/ irefrlayer(inltot), ihypoclayer(inltot), noheadwave,
     &              refraylen(inltot), avhraylen,avvraylen,
     &              hitlay(inltot,3),
     &              irefllayer(inltot),sterr,direrr,refrerr,reflerr
      real avrefrres,avotheres,avreflres,abrefrres,
     & abotheres,abreflres,stnazires
      integer nrrefrres,nrotheres,nrreflres
      common /RESIDUALS/ avrefrres,avotheres,avreflres,
     &                   abrefrres,abotheres,abreflres,
     &                   nrrefrres,nrotheres,nrreflres,
     &                   stnazires(ist,8)
      integer ibackups
      real davar1,xmsqrs1
      common/A3/ ibackups,davar1,xmsqrs1
      real spread
      common/QUALITY/ spread
c
c     Files:
c
      common/files/ modelfilename,stationfilename,seismofilename,
     &             scratchfilename
      character*80 modelfilename,stationfilename,seismofilename,
     &             scratchfilename
      common/infiles/ phasefile,shotfile,topo1file,topo2file,
     &                regnamfile,regkoordfile
      character*80 phasefile,shotfile,topo1file,topo2file,
     &             regnamfile,regkoordfile
      common/outfiles/ velfile,cnvfile,rayfile,outfile,smpfile,stafile,
     &                 drvfile,alefile,dsprfile,rflfile,rfrfile,resfile
      character*80 velfile,cnvfile,rayfile,outfile,smpfile,stafile,
     &             drvfile,alefile,dsprfile,rflfile,rfrfile,resfile
c
c     Various other parameters:
c
      character*80 fm, headerline
      integer      iabort,nitt
      common/FORM/ fm, iabort, nitt, headerline(3)
c
c**************   end of VELESTCOMMON.FOR  *******************************

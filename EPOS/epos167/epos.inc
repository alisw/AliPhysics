c---------------------------------------------------------------------------
c                         dimensions
c---------------------------------------------------------------------------

      parameter (mmry=1)   !memory saving factor

      parameter (mxptl=200000/mmry) !max nr of particles in epos ptl list
      parameter (nmxhep=9990)       !max nr of particles in hep ptl list
      parameter (myptl=1000)        !max nr of droplets in epos ptl list
      parameter (nzeta=60)          !max nr of zeta bins for droplets
      parameter (nflav=6)           !max nr of flavors
      parameter (mxstr=20000/mmry)  !max nr of strings in epos string list 
      parameter (mystr=20000/mmry)    
      parameter (mxtau=4,mxvol=10,mxeps=16)
      parameter (mxtrig=99,mxidco=99)
      parameter (mxpri=200) 
      parameter (mxbins=10000)
      parameter (matau=10,mxcoox=40,mxcooy=10)
      parameter (mxnucl=20)
      parameter (mxhisarg=100)
      parameter (idxD0=0,idxD1=2,idxD=1,nclha=4,nclegy=100) 
      parameter (mamxx=250)
      parameter (mxjerr=10)
           
c---------------------------------------------------------------------------
c                   epos event common block
c---------------------------------------------------------------------------

      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra
      
c     nevt .......... error code. 1=valid event, 0=invalid event
c     bimevt ........ absolute value of impact parameter
c     phievt ........ angle of impact parameter
c     kolevt ........ number of collisions
c     koievt ........ number of inelastic collisions
c     pmxevt ........ reference momentum
c     egyevt ........ pp cm energy (hadron) or string energy (lepton)
c     npjevt ........ number of primary projectile participants
c     ntgevt ........ number of primary target participants
c     npnevt ........ number of primary projectile neutron spectators
c     nppevt ........ number of primary projectile proton spectators
c     ntnevt ........ number of primary target neutron spectators
c     ntpevt ........ number of primary target proton spectators
c     jpnevt ........ number of absolute projectile neutron spectators
c     jppevt ........ number of absolute projectile proton spectators
c     jtnevt ........ number of absolute target neutron spectators
c     jtpevt ........ number of absolute target proton spectators
c     xbjevt ........ bjorken x for dis     
c     qsqevt ........ q**2 for dis  
c     sigtot ........ total cross section
c     nglevt ........ number of collisions acc to  Glauber 
c     zppevt ........ average Z-parton-proj 
c     zptevt ........ average Z-parton-targ

c---------------------------------------------------------------------------
c                   epos particle list common block
c---------------------------------------------------------------------------
      
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)
      common/c1ptl/itsptl(mxptl)
     
c     nptl .......... current particle index (=number of ptls stored)
c     idptl(i) ...... particle id
c     pptl(1,i) ..... x-component of particle momentum 
c     pptl(2,i) ..... y-component of particle momentum 
c     pptl(3,i) ..... z-component of particle momentum 
c     pptl(4,i) ..... particle energy 
c     pptl(5,i) ..... particle mass 
c     iorptl(i) ..... particle number of father (if .le. 0 : no father) 
c     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
c     istptl(i) ..... status: 40 and 41 : Remnant
c                             30 and 31 : Pomeron
c                             20 and 21 : Parton
c                             10 and 11 : Droplet
c                             00 and 01 : Particle
c                            last digit = 0 : last generation
c                            last digit = 1 : not last generation
c     xorptl(1,i) ... x-component of formation point
c     xorptl(2,i) ... y-component of formation point
c     xorptl(3,i) ... z-component of formation point
c     xorptl(4,i) ... formation time
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)
c     ityptl(i)  .... type of particles origin:
c                         10-19: target
c                         20-29: soft Pom
c                         30-39: hard Pom 
c                         40-49: projectile 
c                         50: string, droplet
c     itsptl(i) ..... string type of particles origin (if string)  

      common/c2ptl/iaaptl(mxptl),radptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      common/c4ptl/nptlbd
      common/c6ptl/rinptl(mxptl),vrad,inbxxx
      common/c8ptl/qsqptl(mxptl),zpaptl(mxptl)


c---------------------------------------------------------------------------
c                   hep standard event commonblock.
c---------------------------------------------------------------------------

      double precision phep,vhep

      common/hepevt/nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)
     
c---------------------------------------------------------------------------
c
c         nevhep      -   event number
c         nhep        -   number of entries in the event record
c
c         isthep(i)   -   status code
c         idhep(i)    -   particle id (particle data group standard)
c
c         jmohep(1,i) -   position of mother particle in list
c         jmohep(2,i) -   position of second mother particle in list
c         jdahep(1,i) -   position of first daughter in list
c         jdahep(2,i) -   position of first daughter in list
c
c         phep(1,i)   -   p_x momentum in gev/c
c         phep(2,i)   -   p_y momentum in gev/c
c         phep(3,i)   -   p_z momentum in gev/c
c         phep(4,i)   -   energy in gev
c         phep(5,i)   -   mass in gev/c**2
c
c         vhep(1,i)   -   x position of production vertex in mm 
c         vhep(2,i)   -   y position of production vertex in mm 
c         vhep(3,i)   -   z position of production vertex in mm 
c         vhep(4,i)   -   time of production  in mm/c 
c
c          (note:  1 mm = 10^-12 fm = 5.07 10^-12 1/gev)

c------------------------------------------------------------------------
c  Parameters set in sr aaset and variables to communicate between moduls
c------------------------------------------------------------------------

      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr
      character*80 fnch,fnhi,fndt,fnii,fnid,fnie,fnrj,fnmt
     * ,fngrv,fncp,fnnx,fncs,fndr,fnhy
      common/fname/  fnch, fnhi, fndt, fnii, fnid, fnie, fnrj, fnmt
     * ,fngrv,fncp,fnnx,fncs,fndr,fnhy
      common/nfname/nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhy 
      common/resc2/delvol,deleps,dlzeta,etafac,facnuc,taurea,epscri(3)
      common/resc4/epsdfm
      common/frag1/ndecay,maxres,pud,pudk,pudr,strcut,diqcut
      common/frag2/pdiqua,pdiquak,ioptf,delrex,ptfraqq,ptfrasr,ptfrakd
      common/frag3/aouni,pbreak,pbreakk,ptfra,ptfrair,ptfrak,itflav
      common/frag4/pude,pdiquae,pdiquar,pbreake,pbreakr,ptfrae,ptfrab
      common/frag5/pbreakd,pudd,pdiquad,ptfradr
      common/frag6/qmass(0:6),puds,isospin(0:6)
      real ptq
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      common/hadr10/rstrau(4),rstrad(4),rstras(4),rstrasi
      common/wgtqrk/wgtqqq(4),wgtval,wgtsea,wgtdiq
      double precision timeini,timefin
      common/time1/timeini,timefin
      common/ciotst/iotst1,iotst2,iotst3,iotst4
      common/resc1/taumin,deltau,factau,numtau,amsiac,amprif
      common/resc3/dscale,cepara,iceopt,delamf,deuamf
      common/sprio/ispherio,icotabm,icotabr,icocore
      common/urqm1/iurqmd
      common/incon/cutico,dssico
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      common/othe2/ifrade,iframe,idecay,jdecay
      common/othe3/jframe,kframe
      common/jpsif/jpsi,jpsifi,taumx,nsttau,sigj,ijphis,ijtauan
      common/strlt/iopenu,themas
      common/appli/iappl,model
      common/events/nevent,nfull,nfreeze,ninicon
      common/enrgy/egymin,egymax,elab,ecms,ekin
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/lept1/engy,elepti,elepto,angmue,icinpu
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      common/nucl2/bmaxim,bminim,phimax,phimin
      common/wana1/ymximi,imihis,iclhis,iwtime,wtimet,wtimei,wtimea
      common/wana2/isphis,ispall,wtmini,wtstep,iwcent,iana,nbdky
      common/drop4/asuhax(7),asuhay(7)
      common/gribo/grigam,grirsq,gridel,grislo,gricel,sigppi,sigppd
      common/drop3/bag4rt,dezzer,amuseg,taunll,yradmx,facts,factb,factq
      common/drop2/rcoll,ylongmx,nsegsu,nsegce,facecc,yradpp,yradmi
      common/drop7/ptclu,yradpi,yradpx,ioclude,iocluin,ioquen
      common/metr1/iospec,iocova,iopair,iozero,ioflac,iomom
      common/metr2/nadd,iograc,iocite,ioceau,iociau
      common/hadr2/iomodl,idproj,idtarg,wexcit
      common/hadr25/idprojin,idtargin,rexdifi(4),rexndii(4),irdmpr
      common/metr3/iostat,ioinco,ionlat,ioobsv,iosngl,iorejz,iompar
      common/metr4/ioinfl,ioinct,iowidn,epsgc
      common/lept2/nstmax,prob(99),icbac(99,2),icfor(99,2)
      common/lept3/iolept,igampr,idisco
      common/ebin/noebin,engmin,engmax,nrebin,iologe,iologl
      common/cnsta/pi,pii,hquer,prom,piom,ainfin 
      common/versn/iversn,iverso 
      common/accum/imsg,jerr(mxjerr),ntevt,nrevt,naevt,nrstr,nrptl
      common/accum2/nglacc
      common/cptlu/nptlu /cnrclu/nrclu
      common/drop6/tecm,volu
      common/metr5/iterma,iternc,iterpr,iterpl,iozinc,iozevt
      common/metr6/epsr,keepr 
      common/drop5/keu,ked,kes,kec,keb,ket
      double precision seedi,seedj,seedc
      common/cseed/seedi,seedj,seedc  
      common/cjintc/clust(mxtau,mxvol,mxeps)
      common/cjintd/volsum(mxtau),vo2sum(mxtau),nclsum(mxtau)
      common/ciutot/iutotc,iutote
      common/copen/nopen,nopenr
      common/kopen/kchopen,khiopen,kdtopen,kcpopen,klgopen,knxopen
      character*6 xvaria,yvaria
      common/vana1/xvaria,yvaria,normal,xminim,xmaxim,nrbins,hisfac
      common/vana3/iologb,iocnxb
      parameter(mxnody=200)
      common/nodcy/nrnody,nody(mxnody)
      character*20 subpri
      common/prnt2/nrpri,subpri(mxpri),ishpri(mxpri) 
      common/prnt3/ishevt,ixtau,iwseed,ixgeometry
      parameter (mxcnt=20)
      common/vana4/ar(mxbins,5),ary(mxbins,mxcnt),ardy(mxbins,mxcnt)
     &,ionoerr
      common/xpars/xpar1,xpar2,xpar3,xpar4,xpar5,xpar6,xpar7,xpar8
      common/khist/khisto
      common/ctcor/nctcor/ccttim/ncttim
      common/densi/kdensi(matau,nzeta,mxcoox,mxcooy),tauv(matau)
      common/cjinti/iorsce,iorsdf,iorshh
      common/camim/amimfs,amimel
      common/craddf/scr,scs,hacore
      common/ckoll/iokoll
      common/cncnt/ncnt  /cicnt/inicnt /cnemsi/nemsi
      common/ems1/iemspl,iemsct,gfactor,gwidth
      common/chadron/amproj,amtarg,ypjtl,yhaha,pnullx
      common/vana5/xshift,etacut
      double precision rnucl
      common/nucl5/rnucl(mxnucl,2)
      common/nucl4/nrnucl(2),drnucl(2),rnuclo(mxnucl,2),ionudi
      common/sig/xsig(7),xpom(7)
      common/metr7/ktnbod
      common/hadr3/iregge,isopom,ishpom,iscreen,nprmax,inueff,irmdrop
      common/hadr5/sigtot,sigcut,sigela,sloela,sigsd,sigine,sigdif
     &,sigineaa
      common/hadr6/intpol,isigma,iomega,isetcs
      common/hadr4/alppom,slopom,gamhad(4),r2had(4),chad(4),wdiff(4)
     &            ,gamtil,facdif,r2hads(4),gamhads(4),slopoms,isplit
      common/hadr42/gamhadsi(4)
      common/hadr7/alpreg,sloreg,gamreg,r2reg,ptdiff,ptsend,xmaxremn
     &,xmaxdiff
      common/hadr8/alpqua,alppar,alpsea,alpval,alpdiq,alplea(4),alpdif
      common/hadr14/alpndi,alpdi,ptsendi,zrminc,zrmmax,ptsemx,irzptn
      common/hadr15/zbcut,zopinc,zopmax,zipinc,zipmax,zoeinc,zoemax
      common/hadr16/zidinc,zidmax,zodinc,zodmax,zosinc,zosmax
      common/hadr17/edmaxi,epmaxi
      common/hadr9/ammsqq,ammsqd,ammsdd,cumpom,rexndi(4),rexdif(4)
     &            ,reminv,rexres,zbarfl,rexndf
      common/had10/iclpro,icltar,iclegy
      common/had11/iclpro1,iclpro2,icltar1,icltar2,iclegy1,iclegy2
      common/had12/egylow,egyfac
      common/had13/amdrmax,alpdro(3)
      common/had14/alpcoso,alpcose,betcoso,betcose
      common/emsx1/accept,reject
      common/xems1/iemspr,iemspm,iemspx,iemsrx,iemspu,iemsi2,iemspbx
      common/xems2/iemsse,iemsi1,iemsb,iemsbg,ioems,iemsdr
      common/xspatim/ispacetime
      character*80 hisarg
      common/chisarg/ihisarg,hisarg(2,mxhisarg)
      common /psar10/difnuc(mamxx),radnuc(mamxx)
      parameter (mxbarray=100)
      common/cbarray/barray(mxbarray),nbarray
      common/nxsair/airznxs(3),airanxs(3),airwnxs(3)
     &             ,airavznxs,airavanxs
      common/mod2incs/qgsincs
      common/mod3incs/gheincs
      common/mod4incs/pytincs
      common/mod5incs/hijincs
      common/mod6incs/sibincs
      common/mod7incs/qgsIIincs
      common/testpom/antot,ansh,ansf,pp4max,pp4ini,andropl,anstrg0
     *,anshf,ansff,antotf
     *,anstrg1,anreso0,anreso1,anghadr,antotre
      common/cdiff/anintdiff,anintsdif,anintine
     *,sigineex,sigdifex,sigsdex
      common/cepszer/epszero,alpff(nclha),betff(2)
      common/cgss/tgss(7,7),wgss(7,7)
      
      common/Dparams/alpDs(  idxD0:idxD, nclegy, nclha,nclha),
     &              alpDps( idxD0:idxD, nclegy, nclha,nclha),
     &              alpDpps(idxD0:idxD, nclegy, nclha,nclha),
     &              betDs(  idxD0:idxD, nclegy, nclha,nclha),
     &              betDps( idxD0:idxD, nclegy, nclha,nclha),
     &              betDpps(idxD0:idxD, nclegy, nclha,nclha),
     &              gamDs(  idxD0:idxD, nclegy, nclha,nclha),
     &              delDs(  idxD0:idxD, nclegy, nclha,nclha)
      common/Dparam/alpD(  idxD0:idxD1, nclha, nclha),
     &              alpDp( idxD0:idxD1, nclha, nclha),
     &              alpDpp(idxD0:idxD1, nclha, nclha),
     &              betD(  idxD0:idxD1, nclha, nclha),
     &              betDp( idxD0:idxD1, nclha, nclha),
     &              betDpp(idxD0:idxD1, nclha, nclha),
     &              gamD(  idxD0:idxD1, nclha, nclha),
     &              delD(  idxD0:idxD1, nclha, nclha),
     &     idxDmin,bmxdif(nclha, nclha),bkmxndif
      double precision alpUni,betUni,betpUni,zpUni,ztUni,betfom
      common/DparUni/alpUni(  idxD0:idxD1,2),
     &               betUni(  idxD0:idxD1,2),
     &               betpUni(idxD0:idxD1,2),zpUni,ztUni,
     &               betfom,alpfom,alpfomi,gamfom

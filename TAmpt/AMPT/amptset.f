      SUBROUTINE AMPTSET(EFRM0,FRAME0,PROJ0,TARG0,IAP0,IZP0,IAT0,IZT0)
c
cgsfs added following line to match C++ call
      double precision EFRM0
      double precision xmp, xmu, alpha, rscut2, cutof2
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      CHARACTER*(*) FRAME0,PROJ0,TARG0
      CHARACTER FRAME*8,PROJ*8,TARG*8
      character*25 amptvn
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1A/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
c     parton coalescence radii in case of string melting:
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2
c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
c     initialization value for hadron cascade:
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy
clin-6/2009:
c      common/phidcy/iphidcy
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss
      common/embed/iembed,pxqembd,pyqembd,xembd,yembd
      common/popcorn/ipop

      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
c****************
      EFRM=EFRM0
      FRAME=FRAME0
      PROJ=PROJ0
      TARG=TARG0
      IAP=IAP0
      IZP=IZP0
      IAT=IAT0
      IZT=IZT0

      if(ipop.eq.1) IHPR2(11)=3

clin-6/2009 ctest off turn on jet triggering:
c      IHPR2(3)=1
c     Trigger Pt of high-pt jets in HIJING:
c      HIPR1(10)=7.
c

      if(isoft.eq.1) then
         amptvn = '1.25 (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.25 (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif

      WRITE(*,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#               Version ',a20,             '     #'/1X,
     &10X,'#                06/25/2009                      #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')

c     an odd number is needed for the random number generator:
      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
c     9/26/03 random number generator for f77 compiler:
      CALL SRAND(NSEED)
c
c.....turn on warning messages in nohup.out when an event is repeated:
      IHPR2(10) = 1
c     string formation time:
      ARPAR1(1) = 0.7
c     smearp is the smearing halfwidth on parton z0, 
c     set to 0 for now to avoid overflow in eta.
c     smearh is the smearing halfwidth on string production point z0.
      smearp=0d0
      IAmax=max(iap,iat)
      smearh=1.2d0*IAmax**0.3333d0/(dble(EFRM)/2/0.938d0)
cgsfs Restored this call which was missing
      CALL HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)
      CALL ARTSET
      CALL INIZPC

      RETURN
      END

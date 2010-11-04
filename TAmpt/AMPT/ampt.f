c.....driver program for A Multi-Phase Transport model
      SUBROUTINE AMPT(FRAME0,BMIN,BMAX)
c
      double precision xmp, xmu, alpha, rscut2, cutof2
      double precision smearp,smearh,dpcoal,drcoal,ecritl
cgsfs added following line to match C++ call
      double precision BMIN, BMAX
      integer K
c     CHARACTER*(*) FRAME0
c     CHARACTER FRAME0*8
      CHARACTER*(*) FRAME0
      CHARACTER FRAME*8
cgsfs  added to match specification in AMPTSET
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

      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
c****************

      FRAME=FRAME0
      imiss=0
cgsfs This line should not be here, but the value needs to be set for ARINI2
cgsfs      K=K+1
      K=1

 100  CALL HIJING(FRAME, BMIN, BMAX)
      IAINT2(1) = NATT             


c     evaluate Npart (from primary NN collisions) for both proj and targ:
      call getnp
c     switch for final parton fragmentation:
      IF (IHPR2(20) .EQ. 0) GOTO 2000
c     In the unlikely case of no interaction (even after loop of 20 in HIJING),
c     still repeat the event to get an interaction 
c     (this may have an additional "trigger" effect):
      if(NATT.eq.0) then
         imiss=imiss+1
         if(imiss.le.20) then
            write(6,*) 'repeated event: natt=0,j,imiss=',j,imiss
            goto 100
         else
            write(6,*) 'missed event: natt=0,j=',j
            goto 2000
         endif
      endif
c.....ART initialization and run
      CALL ARINI
      CALL ARINI2(K)
      CALL ARTAN1
      CALL HJANA3
      CALL ARTMN
      CALL HJANA4
      CALL ARTAN2

 2000 CONTINUE
c
c       CALL ARTOUT(NEVNT)
clin-5/2009 ctest off:
c       call flowh0(NEVNT,2)
c       call flowp(2)
c       call iniflw(NEVNT,2)
c       call frztm(NEVNT,2)
c
      RETURN
      END

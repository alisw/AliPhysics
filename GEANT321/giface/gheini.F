*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:14  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.37  by  S.Giani
*-- Author :
      SUBROUTINE GHEINI
C
C *** INITIALIZATION OF RELEVANT GHEISHA VARIABLES ***
C *** INTERFACE WITH GHEISHA8 ***
C *** NVE 20-MAY-1988 CERN GENEVA ***
C
C CALLED BY : GPGHEI, GHEISH
C ORIGIN : F.CARMINATI
C
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/gccuts.inc"
#include "geant321/gsecti.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcking.inc"
#include "geant321/mxgkgh.inc"
C --- GHEISHA COMMONS ---
#include "geant321/s_kginit.inc"
#include "geant321/s_consts.inc"
#include "geant321/s_event.inc"
#include "geant321/s_prntfl.inc"
#include "geant321/s_blank.inc"
#include "geant321/limits.inc"
C
C --- "NEVENT" CHANGED TO "KEVENT" IN COMMON /CURPAR/ DUE TO CLASH ---
C --- WITH VARIABLE "NEVENT" IN GEANT COMMON ---
C
      PARAMETER (MXGKCU=MXGKGH)
      COMMON /CURPAR/ WEIGHT(10),DDELTN,IFILE,IRUN,NEVT,KEVENT,SHFLAG,
     $                ITHST,ITTOT,ITLST,IFRND,TOFCUT,CMOM(5),CENG(5),
     $                RS,S,ENP(10),NP,NM,NN,NR,NO,NZ,IPA(MXGKCU),
     $                ATNO2,ZNO2
C
      DATA CLIGHT /2.99792458E10/
C
C --- SET GHEISHA I/O UNITS TO THE SAME AS FOR GEANT ---
      INBCD=LIN
      NEWBCD=LOUT
C --- CHECK CONSISTENCY OF PREDEFINED VALUES FOR MXGKGH AND MXGKIN.
      IF(MXGKGH .GT. MXGKIN) THEN
         PRINT 1002, MXGKGH,MXGKIN
 1002    FORMAT(1H0,'*** GHEINI ***    MXGKGH = ',I5,' MUST NOT BE ',
     $              'LARGER THAN MXGKIN = ',I5/
     $          1H ,'                  PROGRAM TERMINATED ABNORMALLY')
         STOP
      ENDIF
      IF(MXGKGH .LT. 100)  PRINT 1003, MXGKGH
 1003 FORMAT(1H0,'*** GHEINI ***    WARNING: MXGKGH = ',I5,' SHOULD ',
     $           'BE LARGER THAN 100'/
     $       1H ,'                  PLEASE CHECK')
C
C --- INITIALISE ALL GHEISHA PRINT FLAGS AS FALSE ---
C --- ACTIVATION IS DONE BY "DEBUG" STEERING CARD ---
      NVEFLG=0
      DO 10 J=1,10
      IF (ISWIT(J) .EQ. 109) NVEFLG=1
 10   CONTINUE
      DO 11 J=1,10
      NPRT(J)=.FALSE.
      IF ((J .EQ. 9) .AND. (NVEFLG .EQ. 1)) NPRT(J)=.TRUE.
 11   CONTINUE
      LPRT=.FALSE.
      DO 12 I=1,MXGKPV
         DO 12 J=1,10
            PV(J,I)=0.
  12  CONTINUE
C
C --- INITIALISE KGINIT ARRAY ---
      DO 20 J=1,50
      KGINIT(J)=0
20    CONTINUE
C
C --- INITIALIZE SOME CUT-OFF PARAMETERS WITH GEANT VALUES ---
      TOFCUT=1.0E+20
      NSIZE=MXEVEN
      K0FLAG=0
      CENG(3)=CUTHAD
      CENG(4)=CUTNEU
C
C --- INITIALIZE PI, 2*PI, PI/2 AND PARTICLE PARAMETERS ---
      PI=ACOS(-1.0)
      TWPI=2.0*PI
      PIBTW=PI/2.0
C *** GAMMA ***
      RMASS(1)=0.0
      RCHARG(1)=0.0
      LNVE=LQ(JPART-1)
      RMASS(1)=Q(LNVE+7)
      RCHARG(1)=Q(LNVE+8)
C *** NEUTRINO ***
      RMASS(2)=0.0
      RCHARG(2)=0.0
      LNVE=LQ(JPART-4)
      RMASS(2)=Q(LNVE+7)
      RCHARG(2)=Q(LNVE+8)
C *** E+ ***
      RMASS(3)=0.000511
      RCHARG(3)=1.0
      LNVE=LQ(JPART-2)
      RMASS(3)=Q(LNVE+7)
      RCHARG(3)=Q(LNVE+8)
C *** E- ***
      RMASS(4)=0.000511
      RCHARG(4)=-1.0
      LNVE=LQ(JPART-3)
      RMASS(4)=Q(LNVE+7)
      RCHARG(4)=Q(LNVE+8)
C *** MU+ ***
      RMASS(5)=0.105658
      RCHARG(5)=1.0
      LNVE=LQ(JPART-5)
      RMASS(5)=Q(LNVE+7)
      RCHARG(5)=Q(LNVE+8)
C *** MU- ***
      RMASS(6)=0.105658
      RCHARG(6)=-1.0
      LNVE=LQ(JPART-6)
      RMASS(6)=Q(LNVE+7)
      RCHARG(6)=Q(LNVE+8)
C *** PI+ ***
      RMASS(7)=0.139569
      RCHARG(7)=1.0
      CT=780.4
      LNVE=LQ(JPART-8)
      RMASS(7)=Q(LNVE+7)
      RCHARG(7)=Q(LNVE+8)
      CT=CLIGHT*Q(LNVE+9)
C *** PI0 ***
      RMASS(8)=0.134964
      RCHARG(8)=0.0
      LNVE=LQ(JPART-7)
      RMASS(8)=Q(LNVE+7)
      RCHARG(8)=Q(LNVE+8)
C *** PI- ***
      RMASS(9)=0.139569
      RCHARG(9)=-1.0
      LNVE=LQ(JPART-9)
      RMASS(9)=Q(LNVE+7)
      RCHARG(9)=Q(LNVE+8)
C *** K+ ***
      RMASS(10)=0.493667
      RCHARG(10)=1.0
      CTKCH=370.9
      LNVE=LQ(JPART-11)
      RMASS(10)=Q(LNVE+7)
      RCHARG(10)=Q(LNVE+8)
      CTKCH=CLIGHT*Q(LNVE+9)
C *** K0 SHORT (==> K0) ***
      RMASS(11)=0.49772
      RCHARG(11)=0.0
      CTK0=2.675
      LNVE=LQ(JPART-16)
      RMASS(11)=Q(LNVE+7)
      RCHARG(11)=Q(LNVE+8)
      CTK0=CLIGHT*Q(LNVE+9)
C *** K0 LONG (==> K0 BAR) ***
      RMASS(12)=-0.49772
      RCHARG(12)=0.0
      LNVE=LQ(JPART-10)
      RMASS(12)=-Q(LNVE+7)
      RCHARG(12)=Q(LNVE+8)
C *** K- ***
      RMASS(13)=0.493667
      RCHARG(13)=-1.0
      LNVE=LQ(JPART-12)
      RMASS(13)=Q(LNVE+7)
      RCHARG(13)=Q(LNVE+8)
C *** P ***
      RMASS(14)=0.938272
      RCHARG(14)=1.0
      LNVE=LQ(JPART-14)
      RMASS(14)=Q(LNVE+7)
      RCHARG(14)=Q(LNVE+8)
C *** P BAR ***
      RMASS(15)=-0.938272
      RCHARG(15)=-1.0
      LNVE=LQ(JPART-15)
      RMASS(15)=-Q(LNVE+7)
      RCHARG(15)=Q(LNVE+8)
C *** N ***
      RMASS(16)=0.939566
      RCHARG(16)=0.0
      LNVE=LQ(JPART-13)
      RMASS(16)=Q(LNVE+7)
      RCHARG(16)=Q(LNVE+8)
C *** N BAR ***
      RMASS(17)=-0.939566
      RCHARG(17)=0.0
      LNVE=LQ(JPART-25)
      RMASS(17)=-Q(LNVE+7)
      RCHARG(17)=Q(LNVE+8)
C *** L0 ***
      RMASS(18)=1.11560
      RCHARG(18)=0.0
      CTL0=7.89
      LNVE=LQ(JPART-18)
      RMASS(18)=Q(LNVE+7)
      RCHARG(18)=Q(LNVE+8)
      CTL0=CLIGHT*Q(LNVE+9)
C *** L0 BAR ***
      RMASS(19)=-1.11560
      RCHARG(19)=0.0
      LNVE=LQ(JPART-26)
      RMASS(19)=-Q(LNVE+7)
      RCHARG(19)=Q(LNVE+8)
C *** S+ ***
      RMASS(20)=1.18937
      RCHARG(20)=1.0
      CTSP=2.40
      LNVE=LQ(JPART-19)
      RMASS(20)=Q(LNVE+7)
      RCHARG(20)=Q(LNVE+8)
      CTSP=CLIGHT*Q(LNVE+9)
C *** S0 ***
      RMASS(21)=1.19246
      RCHARG(21)=0.0
      LNVE=LQ(JPART-20)
      RMASS(21)=Q(LNVE+7)
      RCHARG(21)=Q(LNVE+8)
C *** S- ***
      RMASS(22)=1.19734
      RCHARG(22)=-1.0
      CTSM=4.44
      LNVE=LQ(JPART-21)
      RMASS(22)=Q(LNVE+7)
      RCHARG(22)=Q(LNVE+8)
      CTSM=CLIGHT*Q(LNVE+9)
C *** S+ BAR ***
      RMASS(23)=-1.18937
      RCHARG(23)=-1.0
      LNVE=LQ(JPART-27)
      RMASS(23)=-Q(LNVE+7)
      RCHARG(23)=Q(LNVE+8)
C *** S0 BAR ***
      RMASS(24)=-1.19246
      RCHARG(24)=0.0
      LNVE=LQ(JPART-28)
      RMASS(24)=-Q(LNVE+7)
      RCHARG(24)=Q(LNVE+8)
C *** S- BAR ***
      RMASS(25)=-1.19734
      RCHARG(25)=1.0
      LNVE=LQ(JPART-29)
      RMASS(25)=-Q(LNVE+7)
      RCHARG(25)=Q(LNVE+8)
C *** XI0 ***
      RMASS(26)=1.31490
      RCHARG(26)=0.0
      CTX0=8.69
      LNVE=LQ(JPART-22)
      RMASS(26)=Q(LNVE+7)
      RCHARG(26)=Q(LNVE+8)
      CTX0=CLIGHT*Q(LNVE+9)
C *** XI- ***
      RMASS(27)=1.32132
      RCHARG(27)=-1.0
      CTXM=4.92
      LNVE=LQ(JPART-23)
      RMASS(27)=Q(LNVE+7)
      RCHARG(27)=Q(LNVE+8)
      CTXM=CLIGHT*Q(LNVE+9)
C *** XI0 BAR ***
      RMASS(28)=-1.31490
      RCHARG(28)=0.0
      LNVE=LQ(JPART-30)
      RMASS(28)=-Q(LNVE+7)
      RCHARG(28)=Q(LNVE+8)
C *** XI- BAR ***
      RMASS(29)=-1.32132
      RCHARG(29)=1.0
      LNVE=LQ(JPART-31)
      RMASS(29)=-Q(LNVE+7)
      RCHARG(29)=Q(LNVE+8)
C *** DEUTERON ***
      RMASS(30)=1.875613
      RCHARG(30)=1.0
      LNVE=LQ(JPART-45)
      RMASS(30)=Q(LNVE+7)
      RCHARG(30)=Q(LNVE+8)
C *** TRITON ***
      RMASS(31)=2.8144798
      RCHARG(31)=1.0
      LNVE=LQ(JPART-46)
      RMASS(31)=Q(LNVE+7)
      RCHARG(31)=Q(LNVE+8)
C *** ALPHA ***
      RMASS(32)=3.727417
      RCHARG(32)=2.0
      LNVE=LQ(JPART-47)
      RMASS(32)=Q(LNVE+7)
      RCHARG(32)=Q(LNVE+8)
C *** OMEGA- ***
      RMASS(33)=1.67245
      RCHARG(33)=-1.0
      LNVE=LQ(JPART-24)
      RMASS(33)=Q(LNVE+7)
      RCHARG(33)=Q(LNVE+8)
C *** OMEGA- BAR ***
      RMASS(34)=-1.67245
      RCHARG(34)=1.0
      LNVE=LQ(JPART-32)
      RMASS(34)=-Q(LNVE+7)
      RCHARG(34)=Q(LNVE+8)
C *** NEW PARTICLE (GEANTINO) ***
      RMASS(35)=0.0
      RCHARG(35)=0.0
C
      IF (NPRT(9))
     $ PRINT 1000,(I,RMASS(I),RCHARG(I),I=1,33),
     $            CT,CTKCH,CTK0,CTL0,CTSP,CTSM,CTX0,CTXM
 1000 FORMAT(' *GHEINI* === GHEISHA PARTICLE PROPERTIES ==='/
     $ '0INDEX',5X,'MASS (GEV)',5X,'CHARGE'/1H /
     $ 33(1H ,1X,I3,5X,F11.6,6X,F5.2/),
     $ '0PI +-  CT = ',G12.5,' K  +-  CT = ',G12.5/
     $ ' K0     CT = ',G12.5,' L0     CT = ',G12.5/
     $ ' S+     CT = ',G12.5,' S-     CT = ',G12.5/
     $ ' X0     CT = ',G12.5,' X-     CT = ',G12.5)
C
      MP=RMASS(14)
      MPI=RMASS(7)
      MMU=RMASS(5)
      MEL=RMASS(3)
      MKCH=RMASS(10)
      MK0=RMASS(11)
      SMP=MP**2
      SMPI=MPI**2
      SMU=MMU**2
      ML0=RMASS(18)
      MSP=RMASS(20)
      MS0=RMASS(21)
      MSM=RMASS(22)
      MX0=RMASS(26)
      MXM=RMASS(27)
C
C --- LOAD LIMITS FOR INTRINSIC FUNCTION ARGUMENTS ---
      EXPXL = - 82.0
      EXPXU =   82.0
C
      IF (NPRT(9)) PRINT 1001,EXPXL,EXPXU
 1001 FORMAT('0*GHEINI* === INTRINSIC FUNCTION BOUNDARIES ==='/
     $ ' EXPXL,EXPXU = ',2(G12.5,1X))
C
  90  IFINIT(4)=1
C
      END

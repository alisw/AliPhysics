*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.41  by  S.Giani
*-- Author :
*$ CREATE BAMJEV.FOR
*COPY BAMJEV
*
*=== bamjev ===========================================================*
*
      SUBROUTINE BAMJEV ( IHAD, KFA1, KFA2, KFA3, KFA4, AE0, IOPT )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*  Bamjet90: slight revision by A. Ferrari                             *
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
*     Ihad   = number of final hadrons and hadron resonances           *
*     Ae0    = initial energy in GeV                                   *
*     Kfai   = initial quark flavours (u=1,d=2,s=3,c=4,ubar=7,dbar=8,  *
*                                      sbar=9, cbar=10)                *
*     Iopt   = 1,2,3,4,5 means:                                        *
*              1: single (anti)quark jet,   (kfa1)                     *
*              2: single (anti)diquark jet, (kfa1-kfa2)                *
*              3: complete quark antiquark twojet event, (kfa1,kfa2)   *
*              4: complete (anti)quark-(anti)diquark two jet event,    *
*                 (kfa1,kfa2-kfa3)                                     *
*              5: complete diquark-(anti)diquark two jet event,        *
*                 (kfa1-kfa2,kfa3-kfa4)                                *
*     Common/finpar/ contains the momenta,energies and quantum numbers *
*     of the created hadrons                                           *
*     Iv     = actual vertex,iv=1,4,5,6,9,10 are meson verteces        *
*              iv=2,3,7,8 are baryon verteces                          *
*     La     = 1 means cut-off                                         *
*     Ll     = 0,1 means quark jet, antiquark jet, (diquark jet, anti- *
*              diquark jet)                                            *
*     Common/remain/ contains rest jet energy,momenta and quantumnum-  *
*                    bers                                              *
*----------------------------------------------------------------------*
*
#include "geant321/bamjcm.inc"
#include "geant321/finpar2.inc"
#include "geant321/part.inc"
#include "geant321/inpdat.inc"
 
      COMMON/FKREMA/ RPXR,RPYR,RPZR,RER,KR1R,KR2R
      SAVE NCOU, ITMX
*
      DATA NCOU/0/
      DATA ITMX/0/
      NCOU=NCOU+1
      A1SAVE = A1
      B1SAVE = B1
      B2SAVE = B2
      B3SAVE = B3
      IF ( AE0 .LE. 4.D+00 )THEN
         A1 = 0.88D+00
      ELSE IF ( AE0 .LE. 8.0D+00 ) THEN
         A1 = 0.88D+00
      ELSE IF ( AE0 .LT. 30.D+00 ) THEN
         A1 = 0.88D+00
      ELSE
         A1 = 0.88D+00
      END IF
         FRB12 = 0.5D+00
         B1 = 4.D+00 + 1.D+00 / ( FRB12 * AE0 )**2
         B2 = 4.D+00 + 1.D+00 / ( FRB12 * AE0 )**2
      E00 = 40.D+00
* The following is consistent with B3=6
      B3 = B3 * LOG10 ( E00 ) /
     &   ( LOG10 ( 1.D+00 + ( AE0 / E00 )**2 ) + LOG10 ( E00 ) )
C     IF (NCOU.EQ.4701)LT=1
*or   IF (LT.EQ.1)WRITE(LUNOUT,3399)IOPT,NCOU
 3399 FORMAT('  BAMJET',2I10 )
*or   IF (LT.EQ.1)WRITE(LUNOUT,288)IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT
  288 FORMAT (5I5,2E12.4,' BAMJET,IHAD,KFA1,KFA2,KFA3,KFA4,AE0,IOPT')
   63 CONTINUE
      DO 62  I=1,ITMX
         KFR1(I) = 0
         KFR2(I) = 0
   62 CONTINUE
      ITMX=0
   60 CONTINUE
      IYY=0
      IHAD=0
      IT=0
      E0=AE0*.5D0
      IF(IOPT.EQ.1.OR.IOPT.EQ.2) E0=AE0
      LL=0
      IF(KFA1.GT.6.AND.IOPT.EQ.1) LL=1
      IF(KFA1.LE.6.AND.IOPT.EQ.2) LL=1
      IF(KFA1.GT.6.AND.IOPT.EQ.4) LL=1
      PGX = 0.D0
      PGY = 0.D0
      PGZ = 0.D0
      RPX0 = 0.D0
      RPY0 = 0.D0
* The following 6 initializations might be useless, but they make the
* code much clearer
      KR1R = 0
      KR2R = 0
      KR1L = 0
      KR2L = 0
      RER  = 0.D+00
      REL  = 0.D+00
      DO 10 I=1,KMXJCM
         KFR1(I) = 0
         KFR2(I) = 0
         LA = 0
         IT = IT+1
         J  = IT-1
* The following line seems useless
*        K  = IT+1
   40    CONTINUE
C*****CUT OFF TASK
*  | Abbrch is called to cut the chain
         CALL ABBRCH(IT,LL,LA,LT,E0,PGX,PGY,PGZ,KFR1,KFR2,RE,
     &   KR1R,KR2R,KR1L,KR2L,RPX,RPY,RPZ,RPXR,RPYR,RPZR,RPXL,RPYL,RPZL,
     &   RER,REL,IV,B1,B2,KFA1,KFA2,KFA3,KFA4,IOPT,IYY)
         ITMX=MAX(ITMX,IT)
         IF(LA .EQ. 0) GO TO 20
         IT=IT-1
         IF(IOPT.EQ.3.AND.LL.EQ.0) GO TO 70
         IF(IOPT.EQ.4.AND.KFA1.LE.6.AND.LL.EQ.0) GO TO 70
         IF(IOPT.EQ.4.AND.KFA1.GT.6.AND.LL.EQ.1) GO TO 70
         IF(IOPT.EQ.5.AND.LL.EQ.0) GO TO 70
         GO TO 50
   70    CONTINUE
         IYY = 1
         LL  = 1
         IF(IOPT.EQ.4.AND.KFA1.GT.6) LL = 0
         IAR=IT
         GO TO 30
   20    CONTINUE
C*****CHOICE OF THE VERTEX
         CALL FKVERT(IT,LT,LL,KFA1,E0,IV,RE,KFR1,KFR2,AME,IOPT)
C*****CHOICE OF THE FLAVOUR
         CALL FKFLAV(IT,LT,LL,E0,IV,RE,KFR1,KFR2,ISU,BET,KFA1,KFA2,
     &               KFA3,KFA4,IOPT)
C*****CLASSIFICATION OF THE PARTICLES
         CALL HKLASS(IT,LT,LA,LL,KFR1,KFR2,KR1R,KR2R,KR1L,KR2L,IV,IMPS,
     &         IMVE,IB08,IA08,IB10,IA10,AS,B8,KFA1,KFA2,KFA3,KFA4,IOPT)
         ITMX=MAX(ITMX,IT)
         IF (IT .EQ. 1) RX = E0
         IF (IT .GT. 1) RX = RE(J)
         IF(AMF(IT) .GT. RX) GO TO 63
         IF(AMF(IT) .LE. RX) GO TO 19
         LA = 1
         GO TO 40
   19    CONTINUE
         IHAD = IHAD + 1
*or         IF(LT .EQ. 0) GO TO 31
*or         WRITE(LUNOUT,32)IHAD
*or   31    CONTINUE
C*****CHOICE OF THE ENERGY
         CALL ENERGI(IT,LL,LT,IV,RE,HMA,HE,E0,A1)
C*****CHOICE OF THE MOMENTUM
*  |
*  |  He is the total energy, hma the mass one (input) hpx, hpy, hpz
*  |  the momentum components (output values), hps the transversal
*  |  momentum (output)
*  |
         CALL FKIMPU(HE,HMA,HPS,HPX,HPY,HPZ,LT,LL,B3)
         IF (IT .GT. 1) GO TO 13
         RPX(IT)=RPX0-HPX
         RPY(IT)=RPY0-HPY
         GO TO 14
   13    RPX(IT)=RPX(J)-HPX
         RPY(IT)=RPY(J)-HPY
   14    CONTINUE
         IF (IOPT.EQ.1.AND.LL.EQ.1)HPZ=-HPZ
         IF(IOPT.EQ.2.AND.LL.EQ.1) HPZ=-HPZ
         IF(IOPT.EQ.4.AND.KFA1.GT.6) HPZ=-HPZ
         IF(IOPT.EQ.5) HPZ=-HPZ
         PGX=PGX+HPX
         PGY=PGY+HPY
         PGZ=PGZ+HPZ
         PXF(IT)=HPX
         PYF(IT)=HPY
         PZF(IT)=HPZ
*or         IF (LT .EQ. 0) GO TO 15
*or            WRITE(LUNOUT,16)PGX,PGY,PGZ
*or   16       FORMAT(1H0,12HPGX,PGY,PGZ=,3F8.4)
*or   15    CONTINUE
   30    CONTINUE
 
   10 CONTINUE
*
* we suppose that exiting from loop must be achieved via " go to 50
*
      WRITE (LUNERR,*)' BAMJEV: EXITING FROM LOOP ABNORMALLY!!!! '
      WRITE (LUNOUT,*)' BAMJEV: EXITING FROM LOOP ABNORMALLY!!!! '
   50 CONTINUE
      ITMX=MAX(ITMX,IT)
      IF(IOPT.EQ.1.OR.IOPT.EQ.2) GO TO 51
C*****PUT THE RIGHT AND LEFT JET TOGETHER
      CALL VEREIN(IT,LA,LT,RER,REL,RPXR,RPYR,RPZR,RPXL,RPYL,RPZL,
     &KR1R,KR2R,KR1L,KR2L,IHAD,LL,KFR1,KFR2,IMPS,IMVE,IB08,IA08,
     &IB10,IA10,B3,AS,B8,IAR,KFA1,KFA2,KFA3,KFA4,IOPT)
      IF(LA.EQ.3) GO TO 63
      IF(LA.EQ.2) GO TO 63
   51 CONTINUE
      IF(IOPT.EQ.3.OR.IOPT.EQ.4.OR.IOPT.EQ.5) GO TO 52
      IF(LL.EQ.0) GO TO 52
      RPXR=RPXL
      RPYR=RPYL
      RPZR=RPZL
      RER=REL
      KR1R=KR1L
      KR2R=KR2L
   52 CONTINUE
      IF(LE.EQ.0) GO TO 1000
      WRITE(LUNOUT,92)
   92 FORMAT(2X,'NF,NAME,MASS,IQ,IB,PX,PY,PZ,E')
      DO 91 J=1,IT
         WRITE(LUNOUT,90)NREF(J),ANF(J),AMF(J),ICHF(J),IBARF(J),PXF(J),
     &   PYF(J),PZF(J),HEF(J)
   90    FORMAT(2X,I3,A6,F6.3,2I4,4F8.4)
   91 CONTINUE
 1000 CONTINUE
   64 FORMAT(1H0,38HNUMBER OF EVENTS WITH PREST GT. EREST=,I4,
     */,21HNUMBER OF ALL EVENTS=,I4)
 2000 FORMAT(1H0,'NUMBER OF EVENTS WITH ONLY ONE PARTICLE=',I4)
*or      IF(LT.EQ.0) GO TO 17
*or      WRITE(LUNOUT,18)IHAD
*or   18 FORMAT(1H0,15HMULTIPLIZITAET=,I3)
*or   32 FORMAT(1H0,13HHADRONANZAHL=,I3)
*or   17 CONTINUE
      A1 = A1SAVE
      B1 = B1SAVE
      B2 = B2SAVE
      B3 = B3SAVE
      RETURN
      END

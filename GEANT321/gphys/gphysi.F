*
* $Id$
*
* $Log$
* Revision 1.2  1996/09/30 13:28:58  ravndal
* Medium name length checked
*
* Revision 1.1.1.1  1995/10/24 10:21:31  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 06/10/94  16.31.40  by  S.Ravndal
*-- Author :
      SUBROUTINE GPHYSI
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Initialise material constants for all the physics        *
C.    *       mechanisms used by GEANT                                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gclist.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctime.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
      DIMENSION CUTS(10),UCUT(10),MECA(5,13)
      EQUIVALENCE (CUTS(1),CUTGAM),(MECA(1,1),IPAIR)
      CHARACTER*4 DNAME,KCUT(10)
      CHARACTER*20 CHTITL
      LOGICAL NUCRIN
C.
C.    ------------------------------------------------------------------
C.
C              Write RUN parameters, version numbers and CUTS
C
      WRITE(CHMAIL,10000)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10200)GVERSN,IGDATE,IGTIME
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10300)IDRUN
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10400)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10500)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10600)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
C
C              Get the version number of the original INIT structure
C
      OLDGVE=BIG
*
*        Set NUMOLD to 0 to force recalculation of
*        pointers in the tracking routines
      NUMOLD=0
      IF(JRUNG.NE.0)THEN
         OLDGVE =  Q(JRUNG+21)
         IQ(JRUNG+11)=IGDATE
         IQ(JRUNG+12)=IGTIME
         Q(JRUNG+21)=GVERSN
         Q(JRUNG+22)=ZVERSN
C
         DNAME='INIT'
         WRITE(CHMAIL,10700) DNAME,IQ(JRUNG+11),IQ(JRUNG+12), Q(JRUNG+
     +   21), Q(JRUNG+22)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
         DNAME='KINE'
         WRITE(CHMAIL,10700) DNAME,IQ(JRUNG+13),IQ(JRUNG+14), Q(JRUNG+
     +   23), Q(JRUNG+24)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
         DNAME='HITS'
         WRITE(CHMAIL,10700) DNAME,IQ(JRUNG+15),IQ(JRUNG+16), Q(JRUNG+
     +   25), Q(JRUNG+26)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
         DNAME='DIGI'
         WRITE(CHMAIL,10700) DNAME,IQ(JRUNG+17),IQ(JRUNG+18), Q(JRUNG+
     +   27), Q(JRUNG+28)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
         IF(NRNDM(1).EQ.0.AND.NRNDM(2).EQ.0) THEN
*
*             The random number sequence has not been explicitely
*             initialised via a data card. See whether we can initialise
*             it with the words 19/20 of the JRUNG data structure.
            IF(IQ(JRUNG+19).NE.0.OR.IQ(JRUNG+20).NE.0) THEN
               NRNDM(1) = IQ(JRUNG+19)
               NRNDM(2) = IQ(JRUNG+20)
               CALL GRNDMQ(NRNDM(1), NRNDM(2), 0, 'S')
            ENDIF
         ENDIF
         CALL GRNDMQ(IQ(JRUNG+19), IQ(JRUNG+20), 0, 'G')
         WRITE(CHMAIL,10900) IQ(JRUNG+19), IQ(JRUNG+20)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,11000)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
      ENDIF
C
C             Create energy loss and cross-section banks
C
      IF(NEKBIN.LE.0.OR.NEKBIN.GT.199)NEKBIN=90
      IF(EKMIN.GE.EKMAX.OR.EKMIN.LE.0.)THEN
         EKMIN=1.E-5
         EKMAX=1.E+4
      ENDIF
      NEK1=NEKBIN+1
      EKINV=1./LOG10(EKMAX/EKMIN)
      EKBIN(1)=LOG10(EKMIN)
      ELOW(1)=EKMIN
      GEKA=NEKBIN*EKINV
      GEKB=1.-GEKA*EKBIN(1)
      DO 10 I=2,NEK1
         EL=EKBIN(1)+(I-1)/GEKA
         EKBIN(I)=EL
         ELOW(I)=10.**EL
   10 CONTINUE
      ILOW=0
      IF(NMATE.LE.0)GO TO 999
      IF(JMATE.LE.0)GO TO 999
      IF(JTMED.LE.0)GO TO 999
C
      IF(IQ(JTMED-1).LT.40) THEN
         NPUSH=40-IQ(JTMED-1)
         CALL MZPUSH(IXCONS,JTMED,0,NPUSH,'I')
      END IF
      Q(JTMED+31)=ILABS
      Q(JTMED+32)=ISYNC
      Q(JTMED+33)=ISTRA
*
*             If Landau fluctuations activated, cancel delta rays
      KLOS=Q(JTMED+21)
      IF (KLOS .EQ. 0) Q(JTMED+15) = 0.
      IF (KLOS .EQ. 2) THEN
         Q(JTMED+ 8)=9999.
         Q(JTMED+ 9)=9999.
         Q(JTMED+15)=0.
      ENDIF
*
* If Cerenkov generation is on, activate Light absorbtion unless
* explicitely switched off by the user
*
      KLABS=Q(JTMED+31)
      IF(ITCKOV.NE.0) THEN
         IF(KLABS.EQ.-1) THEN
            Q(JTMED+31)=1
         ENDIF
      ENDIF
      Q(JTMED+31)=MAX(Q(JTMED+31),0.)
*
*             If BCUTE,BCUTM,DCUTE,DCUTM,PPCUTM not initialized (=BIG)
*             Set them to CUTGAM,CUTGAM,CUTELE,CUTELE respectively
*
      IF(Q(JTMED+ 6).GT.0.9*BIG)Q(JTMED+ 6)=Q(JTMED+1)
      IF(Q(JTMED+ 7).GT.0.9*BIG)Q(JTMED+ 7)=Q(JTMED+1)
      IF(Q(JTMED+ 8).GT.0.9*BIG)Q(JTMED+ 8)=Q(JTMED+2)
      IF(Q(JTMED+ 9).GT.0.9*BIG)Q(JTMED+ 9)=Q(JTMED+2)
      IF(Q(JTMED+10).GT.0.9*BIG)Q(JTMED+10)=0.010
      IF(Q(JTMED+10).LT.4.*EMASS)Q(JTMED+10)=4.*EMASS
C
      DO 20 K=1,10
   20 CALL GEVKEV(Q(JTMED+K),UCUT(K),KCUT(K))
      WRITE(CHMAIL,10800)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11200)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11300) (UCUT(K),KCUT(K),K=1,3)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11400) (UCUT(K),KCUT(K),K=4,5)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11500) (UCUT(K),KCUT(K),K=6,7)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11600) (UCUT(K),KCUT(K),K=8,10)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11700) (Q(JTMED+K),K=11,13)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,11800) (Q(JTMED+K),K=14,16)
      CALL GMAIL(0,0)
      IF(Q(JTMED+18).EQ.3.) THEN
         NUCRIN = .TRUE.
         Q(JTMED+18)=1.
      ELSE
         NUCRIN = .FALSE.
      ENDIF
      WRITE(CHMAIL,11900) (Q(JTMED+K),K=17,19)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,12000) (Q(JTMED+K),K=20,22)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,12100) Q(JTMED+23),Q(JTMED+31),Q(JTMED+32)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,12110) Q(JTMED+33)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      IF(NUCRIN) THEN
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,12800)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,12900)
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)
         CALL GMAIL(0,0)
      ENDIF
*
* ***    Here we clean up the old cross section tables if any
      DO 40  IMA=1,NMATE
         JMA=LQ(JMATE-IMA)
         IF(JMA.NE.0) THEN
            DO 30  J=1,20
               IF(LQ(JMA-J).NE.0.AND.J.NE.5) THEN
                  CALL MZDROP(IXCONS,LQ(JMA-J),'L')
               ENDIF
   30       CONTINUE
         ENDIF
   40 CONTINUE
*
* *** Call initialisation of the phtotelectric effect constants
      CALL GPHINI
      DO 180 ITM=1,NTMED
         JTM=LQ(JTMED-ITM)
         IF(JTM.LE.0) GO TO 180
         NL=10-IQ(JTM-2)
         IF(NL.GT.0)THEN
            CALL MZPUSH(IXCONS,JTM,NL,0,'I')
            JTM=LQ(JTMED-ITM)
         ENDIF
*         IF(IQ(JTM-1).LT.40) THEN
*           NPUSH=40-IQ(JTM-1)
*           CALL MZPUSH(IXCONS,JTM,0,NPUSH,'I')
*           JTM=LQ(JTMED-ITM)
*         ENDIF
         ISVOL = Q(JTM + 7)
         IFIELD = Q(JTM + 8)
         FIELDM = Q(JTM + 9)
         TMAXFD = Q(JTM + 10)
         STEMAX = Q(JTM + 11)
         DEEMAX = Q(JTM + 12)
         EPSIL = Q(JTM + 13)
         STMIN = Q(JTM + 14)
         IF (TMAXFD.LE.0..OR. (IGAUTO.NE.0.AND.TMAXFD.GT.20.)) THEN
            TMAXFD=20.
            Q(JTM+10) = TMAXFD
         ENDIF
         NMAT = Q(JTM+6)
         JMA = LQ(JMATE-NMAT)
         IF(JMA.LE.0)THEN
            WRITE(CHMAIL,12200)NMAT,ITM
            CALL GMAIL(1,1)
            GO TO 180
         ENDIF
C
C=====>       Get material parameters
C
         A=Q(JMA+6)
         Z=Q(JMA+7)
         DENS=Q(JMA+8)
         RADL=Q(JMA+9)
         IF (Z.LT.1.) THEN
            DEEMAX=0.
            STMIN =0.
            JTP=LQ(JTM)
            IF(JTP.EQ.0) THEN
               CALL MZBOOK(IXCONS,JTP,JTM,0,'TCUT',0,0,40,3,0)
               IQ(JTP-5)=ITM
               DO 50 I=1,23
                  Q(JTP+I)=Q(JTMED+I)
   50          CONTINUE
               Q(JTP+31)=Q(JTMED+31)
               Q(JTP+32)=Q(JTMED+32)
               Q(JTP+33)=Q(JTMED+33)
            ELSEIF(IQ(JTP-1).LT.40) THEN
               NPUSH=40-IQ(JTP-1)
               CALL MZPUSH(IXCONS,JTP,0,NPUSH,'I')
               JTP=LQ(JTM)
               Q(JTP+31)=Q(JTMED+31)
               Q(JTP+32)=Q(JTMED+32)
               Q(JTP+33)=Q(JTMED+33)
            ENDIF
C
C=====>     decay and synch. rad. in vacuum
C
            DO 60 I=11,23
               Q(JTP+I)=0.
   60       CONTINUE
            Q(JTP+20) = Q(JTMED+20)
            Q(JTP+31) = 0.
            Q(JTP+32) = Q(JTMED+32)
            Q(JTP+33) =0.
         ENDIF
C
C=====>       Get tracking medium parameters
C
         JTP=JTMED
         IF(LQ(JTM).NE.0)JTP=LQ(JTM)
         IF(JTP.NE.JTMED)THEN
            IF(IQ(JTP-1).LT.40) THEN
               NPUSH=40-IQ(JTP-1)
               CALL MZPUSH(IXCONS,JTP,0,NPUSH,'I')
               JTP=LQ(JTM)
               Q(JTP+31)=Q(JTMED+31)
               Q(JTP+32)=Q(JTMED+32)
               Q(JTP+33)=Q(JTMED+33)
            ENDIF
            KLOS=Q(JTP+21)
            IF (KLOS .EQ. 2) THEN
               Q(JTP+ 8)=9999.
               Q(JTP+ 9)=9999.
               Q(JTP+15)=0.
            ENDIF
*
* If Cerenkov generation is on, activate Light absorbtion unless
* explicitely switched off by the user
*
            KLABS=Q(JTP+31)
            IF(ITCKOV.NE.0) THEN
               IF(KLABS.EQ.-1) THEN
                  Q(JTP+31)=1
               ENDIF
            ENDIF
            Q(JTP+31)=MAX(Q(JTP+31),0.)
            IF(Q(JTP+ 6).GT.0.9*BIG)Q(JTP+ 6)=Q(JTP+1)
            IF(Q(JTP+ 7).GT.0.9*BIG)Q(JTP+ 7)=Q(JTP+1)
            IF(Q(JTP+ 8).GT.0.9*BIG)Q(JTP+ 8)=Q(JTP+2)
            IF(Q(JTP+ 9).GT.0.9*BIG)Q(JTP+ 9)=Q(JTP+2)
            IF(Q(JTP+10).GT.0.9*BIG)Q(JTP+10)=0.010
            IF(Q(JTP+10).LT.4.*EMASS)Q(JTP+10)=4.*EMASS
*
            CALL UHTOC(IQ(JTM+1),4,CHTITL,20)
            LAST=LNBLNK(CHTITL)
            IF(LAST.GT.0) THEN
               IF(CHTITL(LAST:LAST).EQ.'$') LAST=LAST-1
               IF(LAST.LT.20) CHTITL(LAST+1:20)=' '
            ENDIF
*
            DO 70 K=1,10
   70       CALL GEVKEV(Q(JTP+K),UCUT(K),KCUT(K))
            WRITE(CHMAIL,10100)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,12300)ITM,CHTITL
            CALL GMAIL(0,0)
            WRITE(CHMAIL,12400)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11300) (UCUT(K),KCUT(K),K=1,3)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11400) (UCUT(K),KCUT(K),K=4,5)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11500) (UCUT(K),KCUT(K),K=6,7)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11600) (UCUT(K),KCUT(K),K=8,10)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11700) (Q(JTP+K),K=11,13)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11800) (Q(JTP+K),K=14,16)
            CALL GMAIL(0,0)
            IF(Q(JTP+18).EQ.3.) THEN
               NUCRIN = .TRUE.
               Q(JTP+18)=1.
            ELSE
               NUCRIN = .FALSE.
            ENDIF
            WRITE(CHMAIL,11900) (Q(JTP+K),K=17,19)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,12000) (Q(JTP+K),K=20,22)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,12100) Q(JTP+23),Q(JTP+31),Q(JTP+32)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,12110) Q(JTP+33)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10100)
            CALL GMAIL(0,0)
            IF(NUCRIN) THEN
               WRITE(CHMAIL,10100)
               CALL GMAIL(0,0)
               WRITE(CHMAIL,12800)
               CALL GMAIL(0,0)
               WRITE(CHMAIL,12900)
               CALL GMAIL(0,0)
               WRITE(CHMAIL,10100)
               CALL GMAIL(0,0)
            ENDIF
         ENDIF
C
         DO 80 I=1,10
            CUTS(I)=Q(JTP+I)
   80    CONTINUE
         DO 90 I=1,13
            MECA(1,I)=Q(JTP+10+I)
   90    CONTINUE
         ILABS=Q(JTP+10+21)
         ISYNC=Q(JTP+10+22)
         ISTRA=Q(JTP+10+23)
C
         IF(ILOW.EQ.0)THEN
            DO 100 I=1,10
               IF(Q(JTP+I).LT.0.0000099)THEN
                  WRITE(CHMAIL,12500)
                  CALL GMAIL(1,1)
                  ILOW=1
               ENDIF
  100       CONTINUE
         ENDIF
C
C             Check consistency of different tracking media
C             referencing the same material
C
         DO 120 ITM2=ITM+1,NTMED
            JTM2=LQ(JTMED-ITM2)
            IF(JTM2.NE.0)THEN
               NMAT2=Q(JTM2+6)
               IF(NMAT2.EQ.NMAT)THEN
                  JTP2=JTMED
                  IF(LQ(JTM2).NE.0)JTP2=LQ(JTM2)
                  IF(JTP.NE.JTP2)THEN
                     IF(JTP2.NE.JTMED)THEN
                        KLOS=Q(JTP2+21)
                        IF (KLOS .EQ. 2) THEN
                           Q(JTP2+ 8)=9999.
                           Q(JTP2+ 9)=9999.
                           Q(JTP2+15)=0.
                        ENDIF
                        IF(Q(JTP2+ 6).GT.0.9*BIG)Q(JTP2+ 6)=Q(JTP2+1)
                        IF(Q(JTP2+ 7).GT.0.9*BIG)Q(JTP2+ 7)=Q(JTP2+1)
                        IF(Q(JTP2+ 8).GT.0.9*BIG)Q(JTP2+ 8)=Q(JTP2+2)
                        IF(Q(JTP2+ 9).GT.0.9*BIG)Q(JTP2+ 9)=Q(JTP2+2)
                        IF(Q(JTP2+10).GT.0.9*BIG)Q(JTP2+10)=0.010
                        IF(Q(JTP2+10).LT.4.*EMASS)Q(JTP2+10)=4.*EMASS
                     ENDIF
                     DO 110 I=6,10
                        IF(Q(JTP+I).NE.Q(JTP2+I))THEN
                           WRITE(CHMAIL,12600)NMAT
                           CALL GMAIL(1,0)
                           WRITE(CHMAIL,12700)ITM,ITM2
                           CALL GMAIL(0,1)
                           GO TO 120
                        ENDIF
  110                CONTINUE
                  ENDIF
               ENDIF
            ENDIF
  120    CONTINUE
         IF (DEEMAX.LT.0.) THEN
            IF(ISVOL.EQ.0)THEN
               DEEMAX=0.25
               IF(RADL.GT.2.)DEEMAX=0.25-0.2/SQRT(RADL)
            ELSE
               DEEMAX = 0.2/SQRT(RADL)
            ENDIF
         ENDIF
         IF(OLDGVE.LT.3.15.OR.STEMAX.LE.0.) THEN
*
*       Before version 3.15 there was no STEMAX, so we put it to BIG
            STEMAX=BIG
         ENDIF
         Q(JTM+11) = STEMAX
         Q(JTM+12) = DEEMAX
C
*
*       It can happen that several tracking media refer to the
*       same material. In this case we do not fill the cross section
*       tables more than once. But we still fill the banks of the
*       tracking medium.
         IF(LQ(JMA-1).NE.0) GOTO 160
         NPUSH=20-IQ(JMA-2)
         IF(NPUSH.GT.0)THEN
            CALL MZPUSH(IXCONS,JMA,NPUSH,0,'I')
            JTM=LQ(JTMED-ITM)
            JMA=LQ(JMATE-NMAT)
         ENDIF
*
*     Energy loss and cross-section tables
         IF(ISTRA.EQ.0) THEN
         CALL MZBOOK(IXCONS,LBANK,JMA, -1,'MAEL',0,0,2*NEK1,3,0)
         CALL MZBOOK(IXCONS,LBANK,JMA, -2,'MAMU',0,0, NEK1,3,0)
         ELSE
         CALL MZBOOK(IXCONS,LBANK,JMA, -1,'MAEL',0,0,3*NEK1,3,0)
         CALL MZBOOK(IXCONS,LBANK,JMA, -2,'MAMU',0,0,2*NEK1,3,0)
         ENDIF
         CALL MZBOOK(IXCONS,LBANK,JMA, -3,'MAAL',0,0, NEK1,3,0)
         CALL MZBOOK(IXCONS,JPROB,JMA, -4,'MAPR',0,0, 40,3,0)
         CALL MZBOOK(IXCONS,JPHOT,JMA, -6,'MAPH',2,2, NEK1,3,0)
         CALL MZBOOK(IXCONS,JANNI,JMA, -7,'MAAN',0,0, NEK1,3,0)
         CALL MZBOOK(IXCONS,JCOMP,JMA, -8,'MACO',0,0, NEK1,3,0)
         CALL MZBOOK(IXCONS,JBREM,JMA, -9,'MABR',0,0,3*NEK1,3,0)
         CALL MZBOOK(IXCONS,JPAIR,JMA,-10,'MAPA',0,0,2*NEK1,3,0)
         CALL MZBOOK(IXCONS,JDRAY,JMA,-11,'MADR',0,0,3*NEK1,3,0)
*
* *** Special case for heavy materials, photo-fission
         IF(A.GE.230..AND.A.LE.240..AND.IPFIS.NE.0)THEN
            CALL MZBOOK(IXCONS,JPFIS,JMA,-12,'MAPF',0,0,2*NEK1,3,0)
         ENDIF
*
* *** Rayleigh effect
         CALL MZBOOK(IXCONS,JRAYL,JMA,-13,'MARA',0,0,2*NEK1,3,0)
*
* *** Muon nuclear interactions
         IF(IMUNU.EQ.0)THEN
            JMUNU=0
         ELSE
            CALL MZBOOK(IXCONS,JMUNU,JMA,-14,'MAMN',0,0,NEK1,3,0)
         ENDIF
*
* *** stopping range
         CALL MZBOOK(IXCONS,LBANK,JMA,-15,'MASE',0,0,2*NEK1,3,0)
         CALL MZBOOK(IXCONS,LBANK,JMA,-16,'MASM',0,0,2*NEK1,3,0)
*
* *** Special for photeffect
         CALL GPHXSI
*
* *** coefficients for energy loss
         CALL MZBOOK(IXCONS,LBANK,JMA,-17,'MACE',0,0,6*NEK1,3,0)
         CALL MZBOOK(IXCONS,LBANK,JMA,-18,'MACM',0,0,6*NEK1,3,0)
*
* *** auxiliary tables for integration of dE/dx
         CALL GWORK(NEKBIN*4)
*
         DO 130 JWORK=1, NEKBIN*4
            WS(JWORK) = 0.
  130    CONTINUE
*
* *** Straggling for thin layers, if in effect
         IF(ISTRA.GT.0) THEN
            CALL MZBOOK(IXCONS,JTSTRA,JMA,-19,'MAST',2,2,1,3,0)
#if defined(CERNLIB_ASHO)
            IF(ISTRA.EQ.2) THEN
               CALL MZBOOK(IXCONS,JTASHO,JMA,-20,'MASH',0,0,106,3,0)
            ENDIF
#endif
         ENDIF
*
         DO 140 J=1,20
            JB=LQ(JMA-J)
            IF(JB.NE.0)IQ(JB-5)=NMAT
  140    CONTINUE
C
         JPROB=LQ(JMA-4)
         JMIXT=LQ(JMA-5)
         JPFIS=LQ(JMA-12)
*
* *** Fill above tables (energy losses,cross-sections,stopping ranges)
*
         CALL GPROBI
C
         DO 150 IEKBIN=1,NEK1
C
            CALL GDRELA
            CALL GBRELA
            CALL GPRELA
C
            CALL GPHOTI
            CALL GRAYLI
            CALL GANNII
            CALL GCOMPI
            CALL GBRSGA
            CALL GPRSGA
            CALL GDRSGA
            CALL GMUNUI
            CALL GPFISI
  150    CONTINUE
*
*           Stopping ranges
*
         CALL GRANGI
*
*           Energy loss coefficients
*
         CALL GCOEFF
* *** The table for the energy loss in thin gas layers if the tracking
*     media is defined as such
*
         IF(ISTRA.GT.0) THEN
            CALL GSTINI
#if defined(CERNLIB_ASHO)
            IF (ISTRA.EQ.2) THEN
               CALL GIASHO
            ENDIF
#endif
         ENDIF
*
* *** Multiple scattering,energy-loss and mag.field steps
  160    DO 170 J=1,2
            IF(LQ(JTM-J).NE.0) THEN
               CALL MZDROP(IXCONS,LQ(JTM-J),'L')
            ENDIF
  170    CONTINUE
         CALL MZBOOK(IXCONS,LBANK,JTM, -1,'MUEL',0,0,NEK1+2,3,0)
         IQ(LBANK-5)=ITM
         CALL MZBOOK(IXCONS,LBANK,JTM, -2,'MUMU',0,0,NEK1+2,3,0)
         IQ(LBANK-5)=ITM
         CALL GMULOF
C
  180 CONTINUE
*
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10400)
      CALL GMAIL(0,2)
C
10000 FORMAT(
     +'1************************************************************')
10100 FORMAT(
     +' *                                                          *')
10200 FORMAT(
     +' *    G E A N T  Version',F7.4,'      DATE/TIME',I7,'/',
     +   I4,2X,'*')
10300 FORMAT(
     +' *                      R U N  ',I5,10X,'              *')
10400 FORMAT(
     +' ************************************************************')
10500 FORMAT(
     +' *      Data structure   Date   Time    GVERSN    ZVERSN    *')
10600 FORMAT(
     +' *      --------------   ----   ----    ------    ------    *')
10700 FORMAT(' *',11X,A,6X,I7,2X,I4,3X,F7.4,2X,F7.2,5X,'*')
10800 FORMAT(
     +' *----------------------------------------------------------*')
10900 FORMAT(' *     Random number seeds: ',3X,I10,3X,I10,6X,'*')
11000 FORMAT(
     +' *     --------------------                                 *')
11100 FORMAT(
     +' *              Standard TPAR for this run are              *')
11200 FORMAT(
     +' *              ------------------------------              *')
11300 FORMAT(
     +' *  CUTGAM=',F6.2,A4,'  CUTELE=',F6.2,A4,'  CUTNEU=',F6.2,A4,1X,
     + '*')
11400 FORMAT(
     +' *  CUTHAD=',F6.2,A4,'  CUTMUO=',F6.2,A4,20X,'*')
11500 FORMAT(
     +' *  BCUTE =',F6.2,A4,'  BCUTM =',F6.2,A4,20X,'*')
11600 FORMAT(
     +' *  DCUTE =',F6.2,A4,'  DCUTM =',F6.2,A4,'  PPCUTM=',F6.2,A4,1X,
     + '*')
11700 FORMAT(
     +' *  IPAIR =',F10.0,'  ICOMP =',F10.0,'  IPHOT =',F10.0,1X,'*')
11800 FORMAT(
     +' *  IPFIS =',F10.0,'  IDRAY =',F10.0,'  IANNI =',F10.0,1X,'*')
11900 FORMAT(
     +' *  IBREM =',F10.0,'  IHADR =',F10.0,'  IMUNU =',F10.0,1X,'*')
12000 FORMAT(
     +' *  IDCAY =',F10.0,'  ILOSS =',F10.0,'  IMULS =',F10.0,1X,'*')
12100 FORMAT(
     +' *  IRAYL =',F10.0,'  ILABS =',F10.0,'  ISYNC =',F10.0,1X,'*')
12110 FORMAT(
     +' *  ISTRA =',F10.0, 39X,                                  '*')
 
12200 FORMAT(' ***** GPHYSI error, Material Nr=',I3,
     + ' referenced by Medium Nr=',I3)
12300 FORMAT(
     +' *     Special TPAR for TMED',I4,3X,A,5X,'*')
12400 FORMAT(
     +' *     -------------------------                            *')
12500 FORMAT(' ***** GPHYSI error, CUTS must be',
     + ' greater than 10 KeV *****')
12600 FORMAT(' ***** GPHYSI error for material nr ',I4)
12700 FORMAT(7X,'Tracking medium NR',I4,' and',I4,
     +' have different parameters')
12800 FORMAT(
     +' *  IHADR=3 not supported any more. GHEISHA will handle     *')
12900 FORMAT(
     +' *  hadronic interactions for the above tracking medium     *')
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDECAY
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Control routine for generation of particle decays.            *
C.    *                                                                *
C.    *  MODE(I)        I'th decay mode of current particle            *
C.    *  BRATIO(I)      branching ratio for I'th decay mode of         *
C.    *                 current particle.                              *
C.    *                                                                *
C.    *    ==>Called by : GHSTOP,GTHADR,GTNEUT,GTMUON                  *
C.    *       Author    G.Patrick *********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcking.inc"
#include "geant321/gckine.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcphys.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      DIMENSION BAREA(7)
      DIMENSION BETA(4)
      DIMENSION BRATIO(6)
      DIMENSION MODE(6)
      DIMENSION NTYPE(3)
      DIMENSION PCM(4,3)
      DIMENSION XM(3)
      DIMENSION RNDM(1)
      LOGICAL ROTATE
C.
C.    ------------------------------------------------------------------
C.
C             Search for parent particle in decay list.
C
C
      KCASE   = NAMEC(5)
      NGKINE  = 0
      IF(IDCAY.EQ.2) THEN
         DESTEP = DESTEP+GETOT
         ISTOP  = 2
         GO TO 99
      ENDIF
      DMASS = AMASS
      JPA = LQ(JPART-IPART)
      JDK1 = LQ(JPA-1)
      JDK2 = LQ(JPA-2)
      IF (JDK1.LE.0)                               GO TO 90
      IF (JDK2.LE.0)                               GO TO 90
      DO 5 I=1,6
         BRATIO(I)=Q(JDK1+I)
         MODE(I)=IQ(JDK2+I)
   5  CONTINUE
C
C             Generate branching ratio and select decay mode.
C
      NBR      = 1
      BAREA(1) = 0.
      DO 10 I=2,7
         BRADD    = BRATIO(I-1)
         IF (BRADD.EQ.0.) GO TO 20
         NBR      = NBR+1
         BAREA(I) = BAREA(I-1)+BRADD
  10  CONTINUE
C
  20  CALL GRNDM(RNDM,1)
      BRAND    = 100.*RNDM(1)
      IF (BRAND.GE.BAREA(NBR)) GO TO 99
      ID       = IABS((LOCATF(BAREA,NBR,BRAND)))
C
C             Unpack decay mode.
C
      MXX      = MODE(ID)
      NTYPE(1) = MOD(MXX,100)
      NTYPE(2) = MOD(MXX/100,100)
      JP1 = LQ(JPART-NTYPE(1))
      JP2 = LQ(JPART-NTYPE(2))
      XM(1) = Q(JP1+7)
      XM(2) = Q(JP2+7)
      IF (MXX.LT.10000)THEN
C
C             Two body decay.
C
         NGKINE  = 2
         IF (TLIFE.LT.1.E-15) THEN
            XMTOT = XM(1)+XM(2)
            DO 30 I=1,1000
C--  Create Lorentz distributed energy with FWHM HBAR/TLIFE.
C--  (via integral-transformation of Lorentz-distribution)
C--                 (M.Guckes)
              CALL GRNDM(RNDM,1)
              RMASS = DMASS
     1                + 3.291086E-25/TLIFE * TAN(PI*(RNDM(1)-0.5))
              IF (RMASS.GE.XMTOT) GO TO 40
 30         CONTINUE
            WRITE(CHMAIL,1000) IPART, NTYPE(1), NTYPE(2)
            CALL GMAIL(0,0)
            NGKINE=0
            GO TO 99
 40         DMASS = RMASS
         END IF
         CALL GDECA2(DMASS,XM(1),XM(2),PCM)
      ELSE
C
C             Three body decay.
C
         NTYPE(3) = MXX/10000
         NGKINE  = 3
         JP3 = LQ(JPART-NTYPE(3))
         XM(3) = Q(JP3+7)
         IF (TLIFE.LT.1.E-15) THEN
            XMTOT = XM(1)+XM(2)+XM(3)
            DO 31 I=1,1000
C--  Create Lorentz distributed energy with FWHM HBAR/TLIFE.
C--  (via integral-transformation of Lorentz-distribution)
              CALL GRNDM(RNDM,1)
              RMASS = DMASS
     1                + 3.291086E-25/TLIFE * TAN(PI*(RNDM(1)-0.5))
              IF (RMASS.GE.XMTOT) GO TO 41
 31         CONTINUE
            WRITE(CHMAIL,1000) IPART, NTYPE(1), NTYPE(2), NTYPE(3)
            CALL GMAIL(0,0)
            NGKINE=0
            GO TO 99
 41         DMASS = RMASS
         END IF
         CALL GDECA3(DMASS,XM(1),XM(2),XM(3),PCM)
      ENDIF
C
C             LORENTZ boost into LAB system defined along parent vector
C             followed by rotation back into GEANT system.
C
      P0       = VECT(7)
      E0       = SQRT(P0*P0+DMASS*DMASS)
      BETA(1)  = 0.
      BETA(2)  = 0.
      BETA(3)  = -P0/E0
      BETA(4)  = E0/DMASS
      CALL GFANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
C
      DO 60 K=1,NGKINE
         IF (P0.LE.0.) THEN
            DO 59 I = 1,3
   59       GKIN(I,K) = PCM(I,K)
         ELSE
            CALL GLOREN (BETA, PCM(1,K), GKIN(1,K))
         ENDIF
         IF(ROTATE) CALL GDROT  (GKIN(1,K),COSTH,SINTH,COSPH,SINPH)
         GKIN(4,K)=SQRT(GKIN(1,K)**2+GKIN(2,K)**2+GKIN(3,K)**2+XM(K)**2)
         GKIN(5,K)=NTYPE(K)
         TOFD(K)=0.
         GPOS(1,K) = VECT(1)
         GPOS(2,K) = VECT(2)
         GPOS(3,K) = VECT(3)
  60  CONTINUE
      GO TO 99
C
C             No branching ratio defined. Call user routine
C
#if !defined(CERNLIB_USRJMP)
  90  CALL GUDCAY
#endif
#if defined(CERNLIB_USRJMP)
  90  CALL JUMPT0(JUDCAY)
#endif
C
  99  RETURN
 1000 FORMAT(' ***** GDECAY ERROR : Not enough energy available for ',
     +       'decay of resonance',I3,' to',3I3,'; no decay.')
      END

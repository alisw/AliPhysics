*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.47  by  S.Giani
*-- Author :
      SUBROUTINE GFMDIS
#include "geant321/gcflag.inc"
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
#include "geant321/dimpar.inc"
#if !defined(CERNLIB_SINGLE)
#include "geant321/comcont.inc"
#endif
#include "geant321/comcon.inc"
#if !defined(CERNLIB_SINGLE)
#include "geant321/partt.inc"
#endif
#include "geant321/part.inc"
#include "geant321/gfkdis.inc"
      LOGICAL BTEST
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION SII, ZLL, SIE, ZEL, ONE, PFLUKA
#endif
      PARAMETER (ONE=1)
      DIMENSION IGTOFL(49),IFLTOG(39)
      DATA IGTOFL / 0, 0, 0, 0, 0, 0,23,13,14,12, 15,16, 8, 1, 2,19, 0,
     +17,21,22, 20, 34, 36, 38, 9,18, 31, 32, 33, 35, 37, 39, 17*0/
 
      DATA IFLTOG /14,15, 3, 2, 4, 4, 1,13,25, 5, 6,10, 8, 9,11,12,18,
     +26,16,21, 19,20, 7, 7*0, 27, 28, 29, 22, 30, 23, 31, 24, 32/
      IGF=0
*  Neutrons below 20 MeV kinetic energy passed to MICAP
      IF (IPART.EQ.13.AND.GEKIN.LE.0.02) THEN
         IF (IFINIT(7).EQ.0) CALL GMORIN
* Check that the correct cross-section exists. K.L-P 16.11.93
* BTEST checks if the 0th bit hase been set to 1 (see GMORIN)
         IF (BTEST(IQ(JMA),0))
     +  PRINT *,'*** MICAP: Cross sections for NMAT',NMAT,' not known'
         IGF=2
         SIG = SIGMOR(GEKIN*1.E+9,NMAT)
         IF( SIG .GT. 0.0) THEN
           SHADR = ZINTHA/SIG
         ELSE
           SHADR = BIG
         ENDIF
         GO TO 999
      ENDIF
*
*   FLUKA initialization
      IF (IFINIT(5) .EQ. 0) CALL FLINIT
*   Computation of elastic (SIGEL) and inelastic (NIZLNW)
*   cross-section for each element
      IJ = IGTOFL(IPART)
      ZINE  = BIG
      ZELA  = BIG
      IF(IJ.GT.0) THEN
         PFLUKA = SQRT(GEKIN*(GEKIN+2*AM(IJ)))
         JMA = LQ(JMATE-NMAT)
         NCOMP = ABS(Q(JMA+11))
         DENS = Q(JMA+8)
         JMIXT = LQ(JMA-5)
         IF ( NCOMP .LE. 1) THEN
            CALL NIZLNW(IJ,ONE*Z,ONE*A,ONE*GEKIN,PFLUKA, SII,ZLL)
            IF (ZLL.LT.BIG) THEN
               ZINE = ZLL/DENS
            END IF
         ELSE
            ZIN1 = 0.
            DO 10 K=1,NCOMP
               CALL NIZLNW(IJ,ONE*Q(JMIXT+NCOMP+K),
     +         ONE*Q(JMIXT+K),ONE*GEKIN,PFLUKA,SII,ZLL)
               IF (ZLL.GT.BIG) THEN
                  ZIN1 = 0.0 + ZIN1
               ELSE
                  ZIN1 = DENS*Q(JMIXT+2*NCOMP+K)/ZLL + ZIN1
               END IF
               CABINX(K) = ZIN1
   10       CONTINUE
            ANXNOR = ZIN1
            IF (ZIN1.GT.0.0) THEN
               ZINE = 1./ZIN1
            END IF
         END IF
         IF ( NCOMP .LE. 1) THEN
            CALL SIGEL (IJ,ONE*A,ONE*GEKIN,PFLUKA, SIE,ZEL)
            IF (ZEL.LT.BIG) THEN
               ZELA = ZEL/DENS
            END IF
 
         ELSE
            ZEL1 = 0.
            DO 20 I=1,NCOMP
               CALL SIGEL (IJ,ONE*Q(JMIXT+I),ONE*GEKIN, PFLUKA,
     +         SIE,ZEL)
               IF (ZEL.LT.BIG) THEN
                  ZEL1 = DENS*Q(JMIXT+2*NCOMP+I)/ZEL + ZEL1
               END IF
               CABELX(I) = ZEL1
   20       CONTINUE
            ELXNOR = ZEL1
            IF (ZEL1.GT.0.0) THEN
               ZELA = 1./ZEL1
            END IF
         END IF
      ENDIF
      IF (ZINE.EQ.BIG) THEN
         SINE = 0.0
      ELSE
         SINE = 1./ZINE
      END IF
      IF (ZELA.EQ.BIG) THEN
         SELA = 0.0
      ELSE
         SELA = 1./ZELA
      END IF
      FSIG = SINE + SELA
      IF (FSIG .LE. 0) THEN
         SHADR = BIG
      ELSE
         SHADR = ZINTHA/FSIG
      END IF
  999 CONTINUE
      END

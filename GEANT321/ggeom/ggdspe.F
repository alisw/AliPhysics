*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GGDSPE (JVO, NPAR, PAR, NL, NDIV, ORIG, STEP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GGDSPE (JVO,NPAR,PAR,NL*,NDIV*,ORIG*,STEP*)           *
C.    *                                                                *
C.    *   Computes the actual division parameters of the mother volume *
C.    *    at address JVO with actual parameters in NPAR, PAR          *
C.    *   Returns the division parameters NDIV, ORIG, STEP and the     *
C.    *    number of different division cells NL                       *
C.    *                                                                *
C.    *   Called by : GGDVLP                                           *
C.    *   Author    : S.Banerjee                                       *
C.    *               (Original algorithms of A.McPherson)             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
C.
      PARAMETER (NPAMAX=50)
C.
      DIMENSION PAR(*)
C.
      REAL      PARM(NPAMAX)
      INTEGER   IDTYP(12)
      SAVE IDTYP
      DATA      IDTYP / 1, 1, 1, 1, 2, 2, 3, 3, 4, 1, 6, 6/
C.
C.    ------------------------------------------------------------------
*
* *** Get the division parameters from JVOLUM structure and check if
*         these have to be developed
*
      JIN   = LQ(JVO-1)
      IAXIS = Q(JIN+1)
      NDIV  = Q(JIN+3)
      ORIG  = Q(JIN+4)
      STEP  = Q(JIN+5)
*
      ISHM  = Q(JVO+2)
      IDTY  = IDTYP(ISHM)
      NPARM = Q(JVO+5)
      IF (NPARM.GT.0) THEN
         CALL UCOPY (Q(JVO+7), PARM, NPARM)
      ELSE
         NPARM = NPAR
         CALL VFILL (PARM, NPARM, -1.0)
      ENDIF
*
* *** Find the actual division parameters
*
      IF (IDTY.EQ.1) THEN
*
* BOX, TRD1, TRD2, TRAP, PARA
*
         IF (ISHM.EQ.4) THEN
            IPAR = 1
         ELSE IF (ISHM.EQ.10) THEN
            IPAR = IAXIS
         ELSE
            IPAR = IAXIS + ISHM - 1
         ENDIF
         IF (PARM(IPAR).LT.0.0) ORIG = -PAR(IPAR)
         IF (STEP.LE.0.0) THEN
            STEP = (PAR(IPAR) - ORIG) / NDIV
         ELSE IF (NDIV.LE.0) THEN
            NDIV = (PAR(IPAR) - ORIG + 0.001) / STEP
         ENDIF
         IF (PARM(IPAR).LT.0.0) ORIG = -0.5 * STEP * NDIV
         IF (ISHM.EQ.1.OR.ISHM.EQ.10.OR.(ISHM.EQ.2.AND.IAXIS.EQ.2)) THEN
            NL = 1
         ELSE
            NL = NDIV
         ENDIF
*
      ELSE IF (IDTY.EQ.4) THEN
*
* SPHE
*
         IF (IAXIS.EQ.1.OR.IAXIS.EQ.2) THEN
            IAX1 = 2*IAXIS - 1
            IAX2 = IAX1 + 1
            IF (PARM(IAX1).LT.0.0) ORIG = PAR(IAX1)
            IF (STEP.LE.0.0) THEN
               STEP = (PAR(IAX2) - ORIG) / NDIV
            ELSE IF (NDIV.LE.0) THEN
               NDIV = (PAR(IAX2) - ORIG + 0.001) / STEP
            ENDIF
            IF (PARM(IAX1).LT.0.0) ORIG = 0.5*(ORIG+PAR(IAX2)-STEP*NDIV)
            NL = NDIV
         ELSE
            IF (STEP.LE.0.0.OR.NDIV.LE.0) THEN
               DP = PAR(6) - PAR(5)
               IF (DP.LT.0.0) DP = DP + 360.0
               IF (ORIG.LT.PAR(5)) ORIG = ORIG + 360.0
               IF (ORIG-PAR(5).GT.DP) THEN
                  PMIN = PAR(5)
                  GO TO 910
               ENDIF
               DP = PAR(6) - ORIG
               IF (DP.LT.0.0) DP = DP + 360.0
               IF (NDIV.LE.0) THEN
                  NDIV = (DP + 0.001) /STEP
                  IF (NDIV.LE.0) GO TO 920
               ELSE
                  STEP = DP / NDIV
               ENDIF
            ENDIF
            NL = 1
         ENDIF
*
      ELSE IF (IDTY.EQ.2.OR.IDTY.EQ.3.OR.IDTY.EQ.6) THEN
*
* TUBE, TUBS, CONE, CONS, PGON, PCON
*
         IF (IAXIS.EQ.3) THEN
            IF (IDTY.NE.6) THEN
               IF (IDTY.EQ.2) THEN
                  IPAR = 3
               ELSE
                  IPAR = 1
               ENDIF
               IF (PARM(IPAR).LT.0.0) ORIG = -PAR(IPAR)
               IF (STEP.LE.0.0) THEN
                  STEP = (PAR(IPAR) - ORIG) / NDIV
               ELSE IF (NDIV.LE.0) THEN
                  NDIV = (PAR(IPAR) - ORIG + 0.001) / STEP
               ENDIF
               IF (PARM(IPAR).LT.0.0) ORIG = -0.5 * STEP * NDIV
               IF (IDTY.EQ.2) THEN
                  NL = 1
               ELSE
                  NL = NDIV
               ENDIF
            ELSE
               IF (ISHM.EQ.11) THEN
                  NZ = PAR(4)
               ELSE
                  NZ = PAR(3)
               ENDIF
               IF (NDIV.LE.0.OR.STEP.LE.0.0.OR.NZ.GT.0) GO TO 930
               NL = NDIV
            ENDIF
*
         ELSE IF (IAXIS.EQ.2) THEN
            IF (STEP.LE.0.0.OR.NDIV.LE.0) THEN
               IF (ISHM.EQ.5.OR.ISHM.EQ.7) THEN
                  PMIN = ORIG
                  PMAX = ORIG + 360.0
                  DP   = 360.0
               ELSE IF (ISHM.EQ.6.OR.ISHM.EQ.8) THEN
                  PMIN = PAR(NPAR-1)
                  PMAX = PAR(NPAR)
                  DP   = PMAX - PMIN
               ELSE
                  PMIN = PAR(1)
                  DP   = PAR(2)
                  PMAX = PMIN + DP
               ENDIF
               IF (DP.LT.0.0) DP = DP + 360.0
               IF (ORIG.LT.PMIN) ORIG = ORIG + 360.0
               IF (ORIG-PMIN.GT.DP) GO TO 910
               DP = PMAX - ORIG
               IF (DP.LT.0.0) DP = DP + 360.0
               IF (NDIV.LE.0) THEN
                  NDIV = (DP + 0.001) / STEP
                  IF (NDIV.LE.0) GO TO 920
               ELSE
                  STEP = DP / NDIV
               ENDIF
            ENDIF
            NL = 1
*
         ELSE
            IF (IDTY.NE.6) THEN
               IF (IDTY.EQ.2) THEN
                  IAX1 = 1
                  IAX2 = 2
               ELSE
                  IAX1 = 2
                  IAX2 = 3
               ENDIF
               IF (PARM(IAX1).LT.0.0) ORIG = PAR(IAX1)
               IF (STEP.LE.0.0) THEN
                  STEP = (PAR(IAX2) - ORIG) / NDIV
               ELSE IF (NDIV.LE.0) THEN
                  NDIV = (PAR(IAX2) - ORIG + 0.001) / STEP
               ENDIF
               IF (PARM(IAX1).LT.0.0)
     +            ORIG = ORIG + 0.5 * (PAR(IAX2)-ORIG-STEP*NDIV)
               NL = NDIV
            ELSE
               IF (STEP.LE.0.0.OR.NDIV.LE.0) THEN
                  IF (ISHM.EQ.11) THEN
                     NZ = PAR(4)
                  ELSE
                     NZ = PAR(3)
                  ENDIF
                  GO TO 930
               ENDIF
               NL = NDIV
            ENDIF
         ENDIF
*
      ELSE
         GO TO 900
      ENDIF
      GO TO 999
*
  900 WRITE (CHMAIL, 1001) ISHM, IAXIS
      GO TO 990
*
  910 WRITE (CHMAIL, 1002) ISHM, IAXIS, PMIN, DP, ORIG
      GO TO 990
*
  920 WRITE (CHMAIL, 1003) ISHM, IAXIS, NDIV, DP, STEP
      GO TO 990
*
  930 WRITE (CHMAIL, 1004) ISHM, IAXIS, NDIV, NZ, STEP
*
  990 CALL GMAIL (0, 0)
      IEORUN = 1
*
 1001 FORMAT (' GGDSPE : Invalid call ISHM,IAXIS=',2I5)
 1002 FORMAT (' GGDSPE : Error ISHM,IAXIS,PMIN,DP,ORIG=',2I5,3G12.4)
 1003 FORMAT (' GGDSPE : Error ISHM,IAXIS,NDIV,DP,STEP=',3I5,2G12.4)
 1004 FORMAT (' GGDSPE : Error ISHM,IAXIS,NDIV,NZ,STEP=',4I5,G12.4)
*                                                             END GGDSPE
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GLMOTH (IUDET, IUNUM, NLEV, LVOLS, LINDX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Loads the top part of the Volume tree in LVOLS (IVO's),      *
C.    *   LINDX (IN indices) for a given volume defined through        *
C.    *   its name IUDET and number IUNUM.                             *
C.    *                                                                *
C.    *   The routine stores only upto the last level where JVOLUM     *
C.    *   data structure is developed. If there is no development      *
C.    *   above the current level, it returns NLEV zero.               *
C.    *                                                                *
C.    *   Called by : GDRAW, GFIPAR, GFPARA, GLVOLU                    *
C.    *   Authors   : S.Banerjee                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
C.
      DIMENSION  LVOLS(*), LINDX(*)
      CHARACTER*(*) IUDET
      SAVE MULT
C.
C.    ------------------------------------------------------------------
*
      NLEV = 1
      MULT = IUNUM
      NLDV = 0
      CALL UCTOH (IUDET, LINAM, 4, 4)
      IVOS = IUCOMP (LINAM, IQ(JVOLUM+1), IQ(JVOLUM-1))
      LVOLS(1) = IVOS
*
* *** Scan all volumes to check if it could be mother
*       of the current volume
*
   10 IF (IVOS.GT.1) THEN
         DO 25 IVO = 1, IQ(JVOLUM-1)
            IF (IVO.EQ.IVOS)            GO TO 25
            JVO = LQ(JVOLUM-IVO)
            NIN = Q(JVO+3)
            IF (NIN.EQ.0)               GO TO 25
*
            IF (NIN.LT.0) THEN
               JDIV = LQ(JVO-1)
               IVOT = Q(JDIV+2)
               IF (IVOT.NE.IVOS)        GO TO 25
               IF (MULT.EQ.0) THEN
                  IN   = 1
                  MULT = 1
               ELSE
                  IN   = MULT
               ENDIF
*
            ELSE
               DO 15 IN = 1, NIN
                  JIN  = LQ(JVO-IN)
                  IVOT = Q(JIN+2)
                  IF (IVOT.NE.IVOS)     GO TO 15
                  MULTI = Q(JIN+3)
                  IF (MULT.EQ.0) THEN
                     MULT = MULTI
                  ELSE
                     IF (MULT.NE.MULTI) GO TO 15
                  ENDIF
                  GO TO 20
   15          CONTINUE
               GO TO 25
            ENDIF
*
   20       LINDX(NLEV) = IN
            IF (NLDV.NE.0)              GO TO 30
            NLEV = NLEV + 1
            LVOLS(NLEV) = IVO
            IVOS = IVO
            MULT = 0
            IF (LQ(JVO).NE.0) NLDV = NLEV
            GO TO 10
   25    CONTINUE
      ENDIF
*
* *** Now rearrange the order of the volumes
*
   30 IF (IVOS.EQ.1) LINDX(NLEV) = 1
      NLV  = NLEV/2
      IF (NLV.GT.0) THEN
         DO 35 I = 1, NLV
            J   = NLEV - I + 1
            LID = LINDX(I)
            LVO = LVOLS(I)
            LINDX(I) = LINDX(J)
            LVOLS(I) = LVOLS(J)
            LINDX(J) = LID
            LVOLS(J) = LVO
   35    CONTINUE
      ENDIF
      NLEV = NLEV -1
      IF (LQ(LQ(JVOLUM-LVOLS(1))).EQ.0) NLEV = 0
*                                                             END GLMOTH
  999 END

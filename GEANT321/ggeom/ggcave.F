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
      SUBROUTINE GGCAVE
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   SUBR. GGCAVE                                                 *
C.    *                                                                *
C.    *   Decides whether a given volume is convex or concave.         *
C.    *   For concave volumes, it sets up bit 3 of the status word     *
C.    *   of the corresponding Volume bank                             *
C.    *                                                                *
C.    *   Called by : GGCLOS                                           *
C.    *   Authors   : S.Banerjee                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
*
      IF (IQ(JVOLUM-1).LE.0) GO TO 999
*
* *** Loop over all volumes
*
      DO 20 IVO = 1, IQ(JVOLUM-1)
         JVO = LQ(JVOLUM-IVO)
         ISH = Q(JVO+2)
*
*  **    Some volumes are always convex
*
         IF (ISH.EQ.13)GO TO 20
         IF ((ISH.GT.0.AND.ISH.LE.4).OR.ISH.EQ.10) GO TO 20
*
         IF (ISH.EQ.5.OR.ISH.EQ.6.OR.ISH.EQ.9) THEN
*
*  **       For TUBE/TUBS/SPHR, only when rmin is zero
*  **       For      TUBS/SPHR, D(Phi) should be less than 180.0
*
            IF (Q(JVO+5).GT.0.0) THEN
               IF (Q(JVO+7).EQ.0.0) THEN
                  IF (ISH.EQ.5) THEN
                     GO TO 20
                  ELSE IF (ISH.EQ.6) THEN
                     DPHI = Q(JVO+11) - Q(JVO+10)
                  ELSE
                     DPHI = Q(JVO+12) - Q(JVO+11)
                  ENDIF
                  IF (DPHI.LT.0.0) DPHI = DPHI + 360.0
                  IF (DPHI.LE.180.0.OR.DPHI.GE.360.0) GO TO 20
               ENDIF
            ENDIF
*
         ELSE IF (ISH.EQ.7.OR.ISH.EQ.8) THEN
*
*  **       For CONE/CONS, only when rmin's are zero
*  **       For      CONS, D(Phi) should be less than 180.0
*
            IF (Q(JVO+5).GT.0.0) THEN
               IF (Q(JVO+8).EQ.0.0.AND.Q(JVO+10).EQ.0.0) THEN
                  IF (ISH.EQ.7) THEN
                     GO TO 20
                  ELSE
                     DPHI = Q(JVO+13) - Q(JVO+12)
                     IF (DPHI.LT.0.0) DPHI = DPHI + 360.0
                     IF (DPHI.LE.180.0.OR.DPHI.GE.360.0) GO TO 20
                  ENDIF
               ENDIF
            ENDIF
*
         ELSE IF (ISH.EQ.11.OR.ISH.EQ.12) THEN
*
*  **       For PGON/PCON, only when rmin's are zero, has 2 z-planes
*  **                      and D(phi) less than 180.0 or equal to 360.0
*
            IF (Q(JVO+5).GT.0.0) THEN
               IF (ISH.EQ.11) THEN
                 IPZ  = 10
               ELSE
                 IPZ  = 9
               ENDIF
               NPZ  = Q(JVO+IPZ)
               DPHI = Q(JVO+8)
               IF (DPHI.LT.0.0) DPHI = DPHI + 360.0
               IF (NPZ.EQ.2.AND.Q(JVO+IPZ+2).EQ.0.0.AND.
     +             Q(JVO+IPZ+5).EQ.0.0) THEN
                  IF (DPHI.LE.180.0.OR.DPHI.GE.360.0) GO TO 20
               ENDIF
            ENDIF
*
         ENDIF
C    hyperboloids (ISH=14) are always concave
         IQ(JVO) = IBSET(IQ(JVO),2)
   20 CONTINUE
*                                                             END GGCAVE
  999 END

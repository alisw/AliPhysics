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
      SUBROUTINE GFPARA (NAME, NUMBER, INTEXT, NPAR, NATT, PAR, ATT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Fetch parameters and attributes                             *
C.    *    of the volume NAME-NUMBER.                                  *
C.    *    INTEXT is used to select internal (INTEXT=1)                *
C.    *    or user parameters to be fetched.                           *
C.    *                                                                *
C.    *   Called by : GDPRTR, GDRAW, GDSPEC, GPVOLU                    *
C.    *   Authors   : S.Banerjee, A.McPherson, P.Zanarini              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconsp.inc"
      CHARACTER*4 NAME
      PARAMETER  (NLVMAX=15)
      DIMENSION  PAR(*), ATT(*), LVOLS(NLVMAX), LINDX(NLVMAX)
      LOGICAL BTEST
C.
C.    ------------------------------------------------------------------
C.
      CALL GLOOK (NAME, IQ(JVOLUM+1), NVOLUM, IVO)
      IF (IVO.LE.0) GO TO 991
      JVO  = LQ(JVOLUM-IVO)
      ISH  = Q(JVO+2)
      NPAR = Q(JVO+5)
      NATT = Q(JVO+6)
      JATT = JVO + 6 + NPAR
*
* *** Find a suitable way to access the volume parameters
*
      IF (INTEXT.NE.1.OR..NOT.BTEST(IQ(JVO),1)) THEN
         IF (NPAR.GT.0) THEN
*
*  **       From the JVO structure
*
            JPAR = JVO + 6
         ELSE
*
*  **       From positioning of the volume in the mother
*
            DO 15 IVOM = 1, NVOLUM
               IF (IVO.EQ.IVOM) GO TO 15
               JVOM = LQ(JVOLUM-IVOM)
               NINM = Q(JVOM+3)
               IF (NINM.LE.0)   GO TO 15
               DO 10 IN = 1, NINM
                  JINM = LQ(JVOM-IN)
                  IVOT = Q(JINM+2)
                  IF (IVOT.NE.IVO) GO TO 10
                  NUMR = Q(JINM+3)
                  IF (NUMR.EQ.NUMBER) GO TO 20
   10          CONTINUE
   15       CONTINUE
            GO TO 991
   20       JPAR = JINM + 9
            NPAR = Q(JPAR)
         ENDIF
*
      ELSE
*
*  **    From development structure
*
         CALL GLMOTH (NAME, NUMBER, NLDM, LVOLS, LINDX)
         IF (NLDM.LE.0) GO TO 991
         JVOM = LQ(JVOLUM-LVOLS(NLDM))
         NINM = Q(JVOM+3)
         IF (NINM.LT.0) THEN
            IN   = NUMBER
         ELSE
            DO 25 IN = 1, NINM
               JINM = LQ(JVOM-IN)
               IF (IFIX(Q(JINM+2)).NE.IVO)    GO TO 25
               IF (IFIX(Q(JINM+3)).EQ.NUMBER) GO TO 30
   25       CONTINUE
            GO TO 991
         ENDIF
   30    JPAR = LQ(LQ(JVOLUM-LVOLS(1)))
         IF (NLDM.GT.1) THEN
            DO 35 ILEV = 2, NLDM
               IF (IQ(JPAR+1).EQ.0) THEN
                  JPAR = LQ(JPAR-LINDX(ILEV))
                  IF (JPAR.EQ.0) GO TO 991
               ELSE IF (IQ(JPAR-3).GT.1) THEN
                  JPAR = LQ(JPAR-LINDX(ILEV))
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
   35       CONTINUE
         ENDIF
         IF (NINM.GT.0) THEN
            JPAR = LQ(JPAR-IN)
            IF (JPAR.EQ.0) GO TO 991
         ELSE
            IF (IN.GT.IQ(JPAR+1)) GO TO 991
            IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-IN)
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
         ENDIF
         JPAR = JPAR + 5
         NPAR = IQ(JPAR)
      ENDIF
*
      IF (NPAR.LE.0) GO TO 999
      CALL UCOPY (Q(JPAR+1), PAR, NPAR)
      CALL UCOPY (Q(JATT+1), ATT, NATT)
*
      IF (INTEXT.EQ.1) GO TO 999
      IF (ISH.EQ.28) THEN
*
*  **    NPAR : 30 ---> 12
*
         NPAR = 12
*
      ELSE IF (ISH.EQ.4) THEN
*
*  **    TRAP
*
         NPAR=11
         PH = 0.
         IF (PAR(2).NE.0.) PH = ATAN2(PAR(3),PAR(2))*RADDEG
         TT = SQRT(PAR(2)**2+PAR(3)**2)
         PAR(2) = ATAN(TT)*RADDEG
         IF (PH.LT.0.0) PH = PH+360.0
         PAR(3) = PH
         PAR(7) = ATAN(PAR(7))*RADDEG
         IF (PAR(7).GT.90.0) PAR(7) = PAR(7)-180.0
         PAR(11)= ATAN(PAR(11))*RADDEG
         IF (PAR(11).GT.90.0) PAR(11) = PAR(11)-180.0
*
      ELSE IF (ISH.EQ.10) THEN
*
*  **    PARA
*
         PH = 0.
         IF (PAR(5).NE.0.) PH = ATAN2(PAR(6),PAR(5))*RADDEG
         TT = SQRT(PAR(5)**2+PAR(6)**2)
         PAR(4) = ATAN(PAR(4))*RADDEG
         IF (PAR(4).GT.90.0) PAR(4) = PAR(4)-180.0
         PAR(5) = ATAN(TT)*RADDEG
         IF (PH.LT.0.0) PH = PH+360.0
         PAR(6) = PH
      ENDIF
      GO TO 999
*
  991 NPAR = 0
  999 RETURN
*                                                             END GFPARA
      END

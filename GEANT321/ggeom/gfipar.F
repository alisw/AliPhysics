*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFIPAR (JVO, JIN, IN, NPAR, NATT, PAR, ATT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Routine to fetch internal parameters and attributes for     *
C.    *    the volume from volume address JVO. If it was positioned    *
C.    *    by GSPOSP or declared with negative dimensional parameters, *
C.    *    then JIN and IN must correspond to its positioning inside   *
C.    *    its mother.                                                 *
C.    *                                                                *
C.    *   Called by : GCENT, GFCLIM, GSDVN, GSDVN2, GSDVT, GSDVT2      *
C.    *   Authors   : S.Banerjee, A.McPherson, P.Zanarini              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      PARAMETER  (NLVMAX=15)
      DIMENSION  PAR(*), ATT(*), LVOLS(NLVMAX), LINDX(NLVMAX)
      CHARACTER*4 NAME
      LOGICAL BTEST
C.
C.    ------------------------------------------------------------------
C.
      NPAR = Q(JVO+5)
      NATT = Q(JVO+6)
      JATT = JVO + 6 + NPAR
      IF (.NOT.BTEST(IQ(JVO),1)) THEN
         IF (NPAR.GT.0) THEN
            JPAR = JVO + 6
         ELSE
            IF (JIN.LE.0) GO TO 910
            JPAR = JIN + 9
            NPAR = Q(JPAR)
         ENDIF
      ELSE
         IF (JIN.LE.0) GO TO 910
         JVOM = LQ(JIN+1)
         NIN  = Q(JVOM+3)
         IVO  = Q(JIN+2)
         IF (NIN.GT.0) THEN
            INUM = Q(JIN+3)
         ELSE
            INUM = IN
         ENDIF
         CALL UHTOC(IQ(JVOLUM+IVO),4,NAME,4)
         CALL GLMOTH (NAME, INUM, NLDM, LVOLS, LINDX)
         IF (NLDM.LE.0) GO TO 930
         JPAR = LQ(LQ(JVOLUM-LVOLS(1)))
         IF (NLDM.GT.1) THEN
            DO 10 ILEV = 2, NLDM
               IF (IQ(JPAR+1).EQ.0) THEN
                  JPAR = LQ(JPAR-LINDX(ILEV))
                  IF (JPAR.EQ.0) GO TO 940
               ELSE IF (IQ(JPAR-3).GT.1) THEN
                  JPAR = LQ(JPAR-LINDX(ILEV))
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
   10       CONTINUE
         ENDIF
         IF (NIN.GT.0) THEN
            JPAR = LQ(JPAR-IN)
            IF (JPAR.EQ.0) GO TO 940
         ELSE IF (IQ(JPAR-3).GT.1) THEN
            JPAR = LQ(JPAR-IN)
         ELSE
            JPAR = LQ(JPAR-1)
         ENDIF
         JPAR = JPAR + 5
         NPAR = IQ(JPAR)
      ENDIF
*
      IF (NPAR.LE.0) GO TO 950
      CALL UCOPY (Q(JPAR+1), PAR, NPAR)
      CALL UCOPY (Q(JATT+1), ATT, NATT)
*
      GO TO 999
*
  910 CONTINUE
      WRITE (CHMAIL, 1010) JIN
      CALL GMAIL (0, 0)
 1010 FORMAT (' GFIPAR : Error - JIN = ',I8,' LE 0 where Volume ',
     +        'parameters not present in JVO')
      GO TO 999
*
  930 CONTINUE
      WRITE (CHMAIL, 1030)
      CALL GMAIL (0, 0)
 1030 FORMAT (' GFIPAR : Error - NLDM is zero where development ',
     +        'structure is expected')
      GO TO 999
*
  940 CONTINUE
      WRITE (CHMAIL, 1040)
      CALL GMAIL (0, 0)
 1040 FORMAT (' GFIPAR : Error - JPAR is zero where development ',
     +        'structure is expected')
      GO TO 999
*
  950 CONTINUE
      WRITE (CHMAIL,1050) NPAR
      CALL GMAIL (0, 0)
 1050 FORMAT (' GFIPAR : Error - NPAR = ',I8,' LE zero')
  999 CONTINUE
*                                                             END GFIPAR
      END

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
      SUBROUTINE GMEPOS (JVO, IN, XC, IFL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Updates the common /GCVOLU/ and the structure JGPAR          *
C.    *     for contents positioned.                                   *
C.    *                                                                *
C.    *   For IFL nonzero, it also checks if the point XC is inside    *
C.    *     the content. It returns IN = 0, if the point is outside.   *
C.    *     Otherwise, it transforms XC in the local system.           *
C.    *                                                                *
C.    *   Input : JVO, IN, XC, IFL                                     *
C.    *   Output : IN, XC                                              *
C.    *                                                                *
C.    *   Called by : GDRAW, GMEDIA                                    *
C.    *   Authors   : S.Banerjee, R.Brun, F.Bruyant, A.McPherson       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
      DIMENSION  XC(*)
      REAL       XT(3)
C.
C.    ------------------------------------------------------------------
C.
      JIN  = LQ(JVO-IN)
      IVOT = Q(JIN+2)
      JVOT = LQ(JVOLUM-IVOT)
      IF (NLEVEL.GE.NLDEV(NLEVEL)) THEN
*       (case with JVOLUM structure locally developed)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         DO 10 ILEV = NLDEV(NLEVEL), NLEVEL
            IF (IQ(JPAR+1).EQ.0) THEN
               IF (ILEV.EQ.NLEVEL) THEN
                  JPAR = LQ(JPAR-IN)
               ELSE
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
               ENDIF
               IF (JPAR.EQ.0) GO TO 20
            ELSE IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
   10    CONTINUE
         JPAR = JPAR +5
         NPAR = IQ(JPAR)
         GO TO 30
      ENDIF
*      (normal case)
   20 NPAR = Q(JVOT+5)
      IF (NPAR.EQ.0) THEN
         JPAR = JIN +9
         NPAR = Q(JPAR)
      ELSE
         JPAR = JVOT +6
      ENDIF
*
   30 IROTT  = Q(JIN+4)
      IF (IFL.NE.0) THEN
C*****  Code Expanded From Routine:  GITRAN
C.
C.    ------------------------------------------------------------------
C.
         IF (IROTT .EQ. 0) THEN
            XT(1) = XC(1) - Q(5+JIN)
            XT(2) = XC(2) - Q(6+JIN)
            XT(3) = XC(3) - Q(7+JIN)
*
         ELSE
            XL1 = XC(1) - Q(5+JIN)
            XL2 = XC(2) - Q(6+JIN)
            XL3 = XC(3) - Q(7+JIN)
            JR = LQ(JROTM-IROTT)
            XT(1) = XL1*Q(JR+1) + XL2*Q(JR+2) + XL3*Q(JR+3)
            XT(2) = XL1*Q(JR+4) + XL2*Q(JR+5) + XL3*Q(JR+6)
            XT(3) = XL1*Q(JR+7) + XL2*Q(JR+8) + XL3*Q(JR+9)
*
         ENDIF
C*****  End of Code Expanded From Routine:  GITRAN
         CALL GINME (XT, Q(JVOT+2), Q(JPAR+1), IYES)
         IF (IYES.EQ.0) IN = 0
      ENDIF
*
      IF (IN.GT.0) THEN
*
*       Volume found at deeper level
*
         NL1    = NLEVEL
         NLEVEL = NLEVEL +1
         LVOLUM(NLEVEL) = IVOT
         NAMES(NLEVEL)  = IQ(JVOLUM+IVOT)
         NUMBER(NLEVEL) = Q(JIN+3)
         LINDEX(NLEVEL) = IN
         LINMX(NLEVEL)  = Q(JVO+3)
         GONLY(NLEVEL)  = Q(JIN+8)
         IF (LQ(LQ(JVOLUM-IVOT)).EQ.0) THEN
            NLDEV(NLEVEL) = NLDEV(NL1)
         ELSE
            NLDEV(NLEVEL) = NLEVEL
         ENDIF
         CALL GTRMUL (GTRAN(1,NL1), GRMAT(1,NL1), Q(JIN+5), IROTT,
     +                GTRAN(1,NLEVEL), GRMAT(1,NLEVEL) )
         IQ(JGPAR +NLEVEL) = NPAR
         LQ(JGPAR -NLEVEL) = JPAR
*
         IF (IFL.NE.0) THEN
            XC(1) = XT(1)
            XC(2) = XT(2)
            XC(3) = XT(3)
         ENDIF
      ENDIF
*                                                             END GMEPOS
      END

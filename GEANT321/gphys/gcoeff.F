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
      SUBROUTINE GCOEFF
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates the coefficients for the energy loss               *
C.    *     interpolation                                              *
C.    *     There are 4 tables : electron,positron,muon,proton         *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author      F.Carminati *********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION CX1,CX2,CX3,CY1,CY2,CY3,CDEN1,CDEN2,CDEN3
      DOUBLE PRECISION ACOEFF,BCOEFF,CCOEFF,XRAT,CCOEF1,CCOEF3
      DOUBLE PRECISION SQEPSM,CFACT
      PARAMETER (EPSMAC=1E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1E-11)
#endif
*
      SQEPSM = MAX(1.,91./NEK1)*10.*SQRT(EPSMAC)
      DO 10 IEKBIN=1,NEK1-2
*
         I1 = IEKBIN
         I2 = I1 + 1
         I3 = I2 + 1
         CY1 = ELOW(I1)
         CY2 = ELOW(I2)
         CY3 = ELOW(I3)
         IECOEF = 3*(IEKBIN-1)
*
* *** Electrons
*
         JRANG = LQ(JMA-15)
         JCOEF = LQ(JMA-17)
*
         CX1 = Q(JRANG+I1)
         CX2 = Q(JRANG+I2)
         CX3 = Q(JRANG+I3)
         IF(CX1.NE.CX2.AND.CX1.NE.CX3.AND.CX2.NE.CX3) THEN
            CDEN1 = 1./((CX1-CX2)*(CX1-CX3))
            CDEN2 = 1./((CX2-CX1)*(CX2-CX3))
            CDEN3 = 1./((CX3-CX1)*(CX3-CX2))
            ACOEFF = CY1*CDEN1+CY2*CDEN2+CY3*CDEN3
            BCOEFF = -(CY1*(CX2+CX3)*CDEN1+CY2*(CX1+CX3)*CDEN2+
     +                 CY3*(CX1+CX2)*CDEN3)
            CCOEFF = CY1*CX2*CX3*CDEN1+CX1*CY2*CX3*CDEN2+
     +               CX1*CX2*CY3*CDEN3
            IF(ACOEFF.EQ.0.) THEN
               XRAT = 0.
            ELSEIF(BCOEFF.GT.0.) THEN
               CFACT  = SQRT(ABS(ACOEFF))
               CCOEF1 = SQRT(ABS(CCOEFF-CY1))*CFACT
               CCOEF3 = SQRT(ABS(CCOEFF-CY3))*CFACT
               XRAT   = MAX(CCOEF1,CCOEF3)/BCOEFF
            ELSE
               XRAT=1.
            ENDIF
            IF(XRAT.LE.SQEPSM) THEN
               Q(JCOEF+IECOEF+1) = 0.
               Q(JCOEF+IECOEF+2) = BCOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF
            ELSE
               Q(JCOEF+IECOEF+1) = ACOEFF
               Q(JCOEF+IECOEF+2) = 0.5*BCOEFF/ACOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF/ACOEFF
            ENDIF
         ENDIF
*
* *** Positons
*
         JRANG = LQ(JMA-15) + NEK1
         JCOEF = LQ(JMA-17) +3*NEK1
*
         CX1 = Q(JRANG+I1)
         CX2 = Q(JRANG+I2)
         CX3 = Q(JRANG+I3)
         IF(CX1.NE.CX2.AND.CX1.NE.CX3.AND.CX2.NE.CX3) THEN
            CDEN1 = 1./((CX1-CX2)*(CX1-CX3))
            CDEN2 = 1./((CX2-CX1)*(CX2-CX3))
            CDEN3 = 1./((CX3-CX1)*(CX3-CX2))
            ACOEFF = CY1*CDEN1+CY2*CDEN2+CY3*CDEN3
            BCOEFF = -(CY1*(CX2+CX3)*CDEN1+CY2*(CX1+CX3)*CDEN2+
     +                 CY3*(CX1+CX2)*CDEN3)
            CCOEFF = CY1*CX2*CX3*CDEN1+CX1*CY2*CX3*CDEN2+
     +               CX1*CX2*CY3*CDEN3
            IF(ACOEFF.EQ.0.) THEN
               XRAT = 0.
            ELSEIF(BCOEFF.GT.0.) THEN
               CFACT  = SQRT(ABS(ACOEFF))
               CCOEF1 = SQRT(ABS(CCOEFF-CY1))*CFACT
               CCOEF3 = SQRT(ABS(CCOEFF-CY3))*CFACT
               XRAT   = MAX(CCOEF1,CCOEF3)/BCOEFF
            ELSE
               XRAT=1.
            ENDIF
            IF(XRAT.LE.SQEPSM) THEN
               Q(JCOEF+IECOEF+1) = 0.
               Q(JCOEF+IECOEF+2) = BCOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF
            ELSE
               Q(JCOEF+IECOEF+1) = ACOEFF
               Q(JCOEF+IECOEF+2) = 0.5*BCOEFF/ACOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF/ACOEFF
            ENDIF
         ENDIF
*
* *** Muons
*
         JRANG = LQ(JMA-16)
         JCOEF = LQ(JMA-18)
*
         CX1 = Q(JRANG+I1)
         CX2 = Q(JRANG+I2)
         CX3 = Q(JRANG+I3)
         IF(CX1.NE.CX2.AND.CX1.NE.CX3.AND.CX2.NE.CX3) THEN
            CDEN1 = 1./((CX1-CX2)*(CX1-CX3))
            CDEN2 = 1./((CX2-CX1)*(CX2-CX3))
            CDEN3 = 1./((CX3-CX1)*(CX3-CX2))
            ACOEFF = CY1*CDEN1+CY2*CDEN2+CY3*CDEN3
            BCOEFF = -(CY1*(CX2+CX3)*CDEN1+CY2*(CX1+CX3)*CDEN2+
     +                 CY3*(CX1+CX2)*CDEN3)
            CCOEFF = CY1*CX2*CX3*CDEN1+CX1*CY2*CX3*CDEN2+
     +               CX1*CX2*CY3*CDEN3
            IF(ACOEFF.EQ.0.) THEN
               XRAT = 0.
            ELSEIF(BCOEFF.GT.0.) THEN
               CFACT  = SQRT(ABS(ACOEFF))
               CCOEF1 = SQRT(ABS(CCOEFF-CY1))*CFACT
               CCOEF3 = SQRT(ABS(CCOEFF-CY3))*CFACT
               XRAT   = MAX(CCOEF1,CCOEF3)/BCOEFF
            ELSE
               XRAT=1.
            ENDIF
            IF(XRAT.LE.SQEPSM) THEN
               Q(JCOEF+IECOEF+1) = 0.
               Q(JCOEF+IECOEF+2) = BCOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF
            ELSE
               Q(JCOEF+IECOEF+1) = ACOEFF
               Q(JCOEF+IECOEF+2) = 0.5*BCOEFF/ACOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF/ACOEFF
            ENDIF
         ENDIF
*
* *** Protons
*
         JRANG = LQ(JMA-16) + NEK1
         JCOEF = LQ(JMA-18) +3*NEK1
*
         CX1 = Q(JRANG+I1)
         CX2 = Q(JRANG+I2)
         CX3 = Q(JRANG+I3)
         IF(CX1.NE.CX2.AND.CX1.NE.CX3.AND.CX2.NE.CX3) THEN
            CDEN1 = 1./((CX1-CX2)*(CX1-CX3))
            CDEN2 = 1./((CX2-CX1)*(CX2-CX3))
            CDEN3 = 1./((CX3-CX1)*(CX3-CX2))
            ACOEFF = CY1*CDEN1+CY2*CDEN2+CY3*CDEN3
            BCOEFF = -(CY1*(CX2+CX3)*CDEN1+CY2*(CX1+CX3)*CDEN2+
     +                 CY3*(CX1+CX2)*CDEN3)
            CCOEFF = CY1*CX2*CX3*CDEN1+CX1*CY2*CX3*CDEN2+
     +               CX1*CX2*CY3*CDEN3
            IF(ACOEFF.EQ.0.) THEN
               XRAT = 0.
            ELSEIF(BCOEFF.GT.0.) THEN
               CFACT  = SQRT(ABS(ACOEFF))
               CCOEF1 = SQRT(ABS(CCOEFF-CY1))*CFACT
               CCOEF3 = SQRT(ABS(CCOEFF-CY3))*CFACT
               XRAT   = MAX(CCOEF1,CCOEF3)/BCOEFF
            ELSE
               XRAT=1.
            ENDIF
            IF(XRAT.LE.SQEPSM) THEN
               Q(JCOEF+IECOEF+1) = 0.
               Q(JCOEF+IECOEF+2) = BCOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF
            ELSE
               Q(JCOEF+IECOEF+1) = ACOEFF
               Q(JCOEF+IECOEF+2) = 0.5*BCOEFF/ACOEFF
               Q(JCOEF+IECOEF+3) = CCOEFF/ACOEFF
            ENDIF
         ENDIF
*
   10 CONTINUE
*
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      REAL FUNCTION GPHAK(BET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates Photoelectron Angular ditribution from K shell      *
C.    *   from L1 shell ditribution ENTRY GPHAL1                       *
C.    *                                                                *
C.    *    ==>CALLED BY : GPHOT                                        *
C.    *       AUTHOR    : J. Chwastowski                               *
C.    *                                                                *
C.    ******************************************************************
C.
      DOUBLE PRECISION SAUT1,SAUT2,SAUT3,DRN
      DOUBLE PRECISION COST,SIN2T,X,ANOR,AMAJOR
      DOUBLE PRECISION GAMA,G1,G2,BETA,B1,ONE,HALF
      PARAMETER (BETMAX=0.999999)
      PARAMETER (ONE=1,HALF=ONE/2)
      DIMENSION RAN(2),BETVAL(7),AMAJOR(7)
      DATA BETVAL / 0.8601,0.901,0.9901,0.99901,0.999901,0.9999901,
     +              0.99999901 /
      DATA AMAJOR / 1.02,1.1,4.3,19.,68.,221.,704. /
      DATA ANOR / 0.0 /
C
C L1 shell entry
C
      ENTRY GPHAL1(BET)
C
      BETA = MIN(BET,BETMAX)
      B1 = ONE-BETA
      G2 = ONE/((ONE-BETA)*(ONE+BETA))
      GAMA = SQRT(G2)
      G1 = GAMA-ONE
      DO 10 I = 1,7
       IF(BETA.LE.BETVAL(I)) THEN
        ANOR = AMAJOR(I)*G2
        ANOR = ONE/ANOR
        GO TO 20
       ENDIF
   10 CONTINUE
   20 SAUT2 = HALF*GAMA*G1*(GAMA-2)
      GPHAK = 0.0
   30 CONTINUE
      CALL GRNDM(RAN,2)
      DRN = 2*RAN(1)
      SAUT1 = (DRN-B1)/(DRN*BETA+B1)
      COST = SAUT1
      X = ONE-BETA*COST
      SIN2T = (ONE-COST)*(ONE+COST)
      SAUT3 = SIN2T/(X**2)*(ONE+SAUT2*X)*ANOR
      IF(RAN(2).GT.SAUT3) GO TO 30
      GPHAK = COST
      GPHAL1 = GPHAK
      END

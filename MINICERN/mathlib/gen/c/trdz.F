*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:59  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TRDZ
#include "legbl.inc"
     1,EP
     2,       RGAM,CDIGAM
     3,EP2,PIE,ZZS2,AA,BB,CC,FAC,SUM,ADD
C*NS 2,CLOGOK,RGAM,CDIGAM
C   TRDZ  EXPECTS GR SET=CLOGOK(2.0*ZZ,NFRIG,2)
      K=REAL(VV+1.0)
      EK=K
      EP=VV+0.5-EK
      ZZS2=Z2/4.0
      KK=K-1
      SUM=0.0
      IF(K)9,9,10
   10 B=0.5
      C=1.0
      FAC=RGAM(VV,B,C)*EXP(VV*GR)
      SUM=FAC
      IF(KK)9,11,12
   12 AA=-VV-2.0
      BB=-VV-0.5
      DO 2 I=1,KK
      BB=BB+1.0
      AA=AA+2.0
      FAC=FAC*(1.0+AA)*AA*ZZS2/(BB*I)
      SUM=SUM+FAC
      IF(CRIT(SUM,FAC,ACC))51,2,2
    2 CONTINUE
      GOTO 11
    9 SUM=0.0
   51 CONTINUE
   11 PIE=PI*EP
      EP2=EP/2.0
      A=EK
      B=0.5-EP
      C=1.0
      BB=0.0
      FAC=-(1.0-PIE**2/3.0)*RGAM(A,B,C)*RGAM(C,BB,-EP)
     1                 / (PI*EXP((2.0*EK-VV)*GR))
      A=EK+0.5
      B=1.0-EP2
      C=EK+1.0+EP2
      ADD=2.*CDIGAM(A)-CDIGAM(B)-CDIGAM(C)-2.0*GR
      GR=FAC*ADD*(1.0+EP2*ADD)
      SUM=SUM+GR
      AA=EK-1.5-EP
      BB=-EP
      CC=EK
      DO 22 I=1,50
      AA=AA+2.0
      BB=BB+1.0
      CC=CC+1.0
      FAC=FAC*(1.0+AA)*AA*ZZS2/(BB*CC)
      ADD=ADD+2.0/(A+1.0)+2.0/A-1.0/B-1.0/C
      A=A+2.0
      B=B+1.0
      C=C+1.0
      GR=FAC*ADD*(1.0+EP2*ADD)
      SUM=SUM+GR
      IF(CRIT(SUM,GR,ACC))52,22,22
   22 CONTINUE
      NCVG=NCVG+4
   52 PP=SUM/PISR
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE TWOPAR.FOR
*COPY TWOPAR
*
*=== twopar ===========================================================*
*
      SUBROUTINE TWOPAR(E1,E2,P1,P2,COD1,COD2,COF1,COF2,SIF1,SIF2,
     1IT1,IT2,UMOO,ECM,P,N,AM1,AM2)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C     ******************************************************
C     QUASI TWO PARTICLE PRODUCTION
C     TWOPAR CALCULATES THE ENERGYS AND THE MOMENTA
C     FOR THE CREATED PARTICLES OR RESONANCES IT1 AND IT2
C     IN THE CM - SYSTEM
C     COD1,COD2,COF1,COF2,SIF1,SIF2 ARE THE ANGLES FOR
C     SPHERICAL COORDINATES
C     ******************************************************
      COMMON /FKGAMR/ REDU,AMO,AMM(15 )
      COMMON/FKABLT/AM(110),GA(110),TAU(110),ICH(110)
     *,IBAR(110),K1(110),K2(110)
   9   CONTINUE
      AMA=AM1
      AMB=AM2
      AMA2=AMA*AMA
      AMB2=AMB*AMB
      S=UMOO*UMOO
      E1=(S-AMB2+AMA2)/(2.D0*UMOO)
      E2=(S-AMA2+AMB2)/(2.D0*UMOO)
      AMTE=MAX(E1*E1-AMA2,ZERZER)
  11  CONTINUE
      AMTE=AMTE
      P1=SQRT(AMTE)
      P2=P1
C     / P2 / = / P1 /  BUT OPPOSITE DIRECTIONS
C     DETERMINATION  OF  THE ANGLES
C     COS(THETA1)=COD1      COS(THETA2)=COD2
C     SIN(PHI1)=SIF1        SIN(PHI2)=SIF2
C     COS(PHI1)=COF1        COS(PHI2)=COF2
C     PHI IS UNIFORMLY DISTRIBUTED IN ( 0,2*PI )
      CALL COSI(COF1,SIF1)
      COF2=-COF1
      SIF2=-SIF1
C     CALCULATION OF THETA1
      CALL TCHOIC(TR,P,P1,ECM,E1,IT1,IT2,N,AM1,AM2)
      COD1=(TR-AMA2-AM(N)*AM(N)+2.D0*ECM*E1)/(2.D0*P*P1+1.D-18)
      IF (ABS(COD1) .GT. 1.D+00) COD1 = SIGN (ONEONE,COD1)
      COD2=-COD1
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE ENERGI.FOR
*COPY ENERGI
*
*=== energi ===========================================================*
*
      SUBROUTINE ENERGI(IT,LL,LT,IV,RE,HMA,HE,E0,A1)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*  Energi89: slight revision by A. Ferrari                             *
*----------------------------------------------------------------------*
*
*      INCLUDE '(BAMJCM)'
#include "geant321/finpar2.inc"
      DIMENSION IV(*),RE(*)
      REAL RNDM(1)
C     SAMPLES THE ENERGY OF THE CREATED HADRON
C*****HE IS THE ENERGY OF THE ACTUAL HADRON
      I=IT
      J=IT-1
      HMA=AMF(I)
      IF(IT .EQ. 1) RX=E0
      IF(IT .GT. 1) RX=RE(J)
      CALL GRNDM(RNDM,1)
      X   = RNDM(1)
      IVY = IV(I)
      GO TO (1,2,2,1,1,1,2,2,1,1),IVY
C     MESON VERTEX
    1 AX=(1.D0-A1)/(3.D0*A1)
      C=1.D0/((1.D0-A1)*(1.D0-HMA/RX)+3.D0*.5D0*A1*(1.D0-HMA/RX)**2)
      ETA=-AX+SQRT(AX**2+X*2.D0/(3.D0*A1*C))
      GO TO 3
C     BARYON VERTEX
    2 CY=1.D0/((1.D0-A1)*(1.D0-HMA/RX)+A1*(1.D0-HMA/RX)**3)
      BY=X/(2.D0*A1*CY)
      BZ=((1.D0-A1)/(3.D0*A1))**3
      BV=BY+SQRT(BY**2+BZ)
      BU=BY-SQRT(BY**2+BZ)
      IF(BU.LT.0.D0) BU1=ABS(BU)
      BU3=BU1**(1.D0/3.D0)
      BV3=BV**(1.D0/3.D0)
      IF(BU.LT.0.D0) BU3=-BU3
      ETA=BU3+BV3
      GO TO 3
    3 CONTINUE
      RE(I)=ETA*RX
      HE=(1.D0-ETA)*RX
      HEF(I)=HE
      RETURN
      END

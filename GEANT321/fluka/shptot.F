*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE SHPTOT.FOR
*COPY SHPTOT
*
*=== shptot ===========================================================*
*
      FUNCTION SHPTOT(IT,PO)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C   LAST CHANGE 25.11 -86 BY PERTTI AARNIO
C
C
C   TOTAL HADRON-PROTON CROSS SECTIONS
C   PLAB.GE.10 GEV
C
C********************************************************************
#include "geant321/paprop.inc"
      REAL RNDM(1)
C
      F1=1.D0
      ITT=IT
      AMIT2=AM(ITT)**2
      UMO2=AMIT2 + AM(1)**2 + 2.D0*AM(1)*(PO+0.5D0*AMIT2/PO)
      UMO=SQRT(UMO2)
C
      A4=0.D0
      A5=0.D0
      A6=0.D0
C
      GO TO
     * (1,1,50,50,50,50,50,2,2,50,50,5,3,3,4,4,8,8,5,8,8,8,50,6,7,
     *  1,1,1,1,1,1,1,1,1,1,1,1,1,1), ITT
C
 1    CONTINUE
      A1=38.4D0
      A2=0.46D0
      A3=125.D0
      IF(ITT.EQ.1) GOTO 100
      A5=84.1D0
      A6=0.43D0
      GOTO 100
C
 2    CONTINUE
      A1=38.5D0
      A2=0.46D0
      A3=125.D0
      A4=15.D0
      IF(ITT.EQ.8) GOTO 100
      A5=77.43D0
      A6=0.40D0
      GOTO 100
C
 3    CONTINUE
      IF(UMO.LT.47.D0) GOTO 31
      F1=0.6667D0
      ITT=1
      GOTO 1
 31   CONTINUE
      A1=24.D0
      A2=0.60D0
      A3=160.D0
*   Very very crude patch for the pion- elastic xsec problem
*     IF(ITT.EQ.13) GOTO 100
      IF(ITT.EQ.13.OR.ITT.EQ.14) GOTO 100
      A5=7.9D0
      A6=0.54D0
      GOTO 100
C
 4    CONTINUE
      IF(UMO.LT.110.D0) GOTO 41
      F1=0.6667D0
      ITT=1
      GOTO 1
 41   CONTINUE
      A1=20.3D0
      A2=0.59D0
      A3=140.D0
      IF(ITT.EQ.15) GOTO 100
      A5=30.13D0
      A6=0.42D0
      GOTO 100
C
 5    CONTINUE
      ITT=15
      CALL GRNDM(RNDM,1)
      IF(RNDM(1).LT.0.5D0) ITT=16
      GOTO 4
C
 6    CONTINUE
C***
C   K-ZERO:  SET EQUAL TO K-/PROTON
C            (SHOULD BE K-/NEUTRON)
C***
      ITT=16
      GOTO 4
C
 7    CONTINUE
C***
C   K-ZERO BAR:  SET EQUAL TO K+/PROTON
C                (SHOULD BE K+/NEUTRON)
C***
      ITT=15
      GOTO 4
C
 8    CONTINUE
C***
C   SIGMA +/-/0  AND  LAMBDA/LAMBDA BAR:  SET EQUAL TO P-P
C***
      ITT=1
      GOTO 1
C
 50   CONTINUE
C***
C   LEPTONS AND PI0
C***
      SHPTOT=1.D-10
      RETURN
C
 100  CONTINUE
C
      SHPTOT=A1+A2*(LOG(UMO2/A3))**2+A4/UMO2+A5*UMO2**A6
      SHPTOT=F1*SHPTOT
      RETURN
      END

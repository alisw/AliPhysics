*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE XLAMB.FOR
*COPY XLAMB
*
*=== xlamb ============================================================*
*
      FUNCTION XLAMB(X,Y,Z)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      COMMON /FKIDGB/ IDGB,IDG
      COMMON /FKGAMR/ REDU,AMO,AMM(15 )
      COMMON/FKDREI/TEST(12)
      COMMON/FKPRUN/ISYS
      YZ=Y-Z
      XLAMB=X*X-2.D0*X*(Y+Z)+YZ*YZ
      XLAM =XLAMB
      IF (IDGB.LE.0) GO TO 11
      IF(XLAMB.GT.1.D-12) GOTO 11
      WRITE(ISYS,12)
      WRITE(ISYS,10) XLAM,X,Y,Z,TEST
      WRITE(ISYS,13)
   12 FORMAT(/,10X,12H XLAMB PRINT)
   13 FORMAT(10X,60(1H*))
   10 FORMAT(4E20.8,5HXLAMB,/,12F10.5)
   11 CONTINUE
      IF(XLAMB.LE.0.D0)XLAMB=ABS(XLAMB)
      XLAMB=SQRT(XLAMB)
      RETURN
      END

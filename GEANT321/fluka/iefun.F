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
*$ CREATE IEFUN.FOR
*COPY IEFUN
*
*=== iefun ============================================================*
*
      INTEGER FUNCTION IEFUN(PL,IRE)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/reac.inc"
#include "geant321/redver.inc"
*
C*****IEFUN CALCULATES A MOMENTUM INDEX
      COMMON / FKRUN    / RUNTES, EFTES
*
      IPLA=IEII(IRE)+1+1
      IPLE=IEII(IRE+1)
      IF (PL.LT.0.D0) GO TO 13
      DO 2 I=IPLA,IPLE
      J=I-IPLA+1
      IF (PL.LE.PLABF(I)) GO TO 3
    2 CONTINUE
      I=IPLE
      IF ( EFTES.GT.40.D0) GO TO 4321
      EFTES=EFTES+1.D0
      WRITE(LUNOUT,10)PL,J
 4321 CONTINUE
      GO TO 4323
   13 CONTINUE
      DO 12 I=IPLA,IPLE
      J=I-IPLA+1
      IF (-PL.LE.UMO(I)) GO TO 3
   12 CONTINUE
      I=IPLE
      IF ( EFTES.GT.40.D0) GO TO 4322
      EFTES=EFTES+1.D0
      WRITE(LUNOUT,10)PL,I
 4322 CONTINUE
    3 CONTINUE
 4323 CONTINUE
      IEFUN=I
      RETURN
   10 FORMAT(14H PLAB OR -ECM=,E10.4,27H IS OUT OF CONSIDERED RANGE
     *,7H IEFUN=,I5)
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE CHANWT.FOR
*COPY CHANWT
*
*=== CHANWT ===========================================================*
*
      SUBROUTINE CHANWT
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
#include "geant321/decayc.inc"
#include "geant321/part.inc"
      COMMON/FKPRUN/ISYS
      DIMENSION HWT(IDMAX9)
C  CHANGE OF WEIGHTS WT FROM ABSOLUT VALUES INTO THE SUM OF WT OF A DEC.
      DO 3  J=1,IDMAX9
    3 HWT(J)=0.D0
      DO 1  I=1,IDMAX8
      IK1=K1(I)
      IK2=K2(I)
      HV=0.D0
      DO 2 J=IK1,IK2
      HV=HV+WT(J)
      HWT(J)=HV
*
*  Modified by A.Ferrari to avoid false error messages
*
*     IF(HWT(J).GT.1.D0) WRITE(ISYS,101) HWT(J),J,I,IK1
      IF(HWT(J) - 1.D0 .GT. 1.D-14) WRITE(ISYS,101) HWT(J),J,I,IK1
    2 CONTINUE
  101 FORMAT(2X,15H ERROR IN HWT =,1F10.5,8H J,I,K1=,3I5)
    1 CONTINUE
      DO 4  J=1,IDMAX9
    4 WT(J)=HWT(J)
      RETURN
      END

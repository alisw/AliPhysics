*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE D501P2(LRP,N,A,B,C,LAMU,PHI,DPHINO,ITER,LFN,MODE,VERS)

************************************************************************
*   LEAMAX, VERSION: 15.03.1993
************************************************************************
*
*   THIS ROUTINE CONTROLS THE PRINTING OF THE PACKAGE LEAMAX.
*
************************************************************************

#include "gen/imp64.inc"
#include "gen/def64.inc"
     +   LAMBDA,LAMU
      LOGICAL LFN,LRP
      CHARACTER VERS*6,TIT(2)*18
      DIMENSION A(*),B(*),C(*),LAMU(*)


      DATA TIT(1),TIT(2) /' ','STANDARD DEVIATION'/

      IF(.NOT.LRP) THEN

       WRITE(6,1030)
       LRP=.TRUE.
      ENDIF

      IF(LFN) THEN
       WRITE(6,1010) 'END:',ITER,PHI,DPHINO
       IF(VERS .EQ. 'DFUNFT') THEN
        WRITE(6,1020) TIT(2)(1:8),TIT(2)(10:18)
       ELSE
        WRITE(6,1020) TIT(1)(1:8),TIT(1)(10:18)
       ENDIF
      ELSE
       WRITE(6,1010) '    ',ITER,PHI,DPHINO
       WRITE(6,1020) TIT(1)(1:8),TIT(1)(10:18)
      ENDIF

      IF(LFN .AND. VERS .EQ. 'DFUNFT') THEN
       WRITE(6,1040) (I,A(I),B(I),LAMU(I),C(I), I=1,N)
      ELSE
       WRITE(6,1050) (I,A(I),B(I),LAMU(I), I=1,N)
      ENDIF

      RETURN

 1010 FORMAT(/6X,A4,' ITERATION',I5,3X,'PHI = ',1PD12.5,6X,
     1       'GNO = ',1PD12.5/)
 1020 FORMAT(12X,'PARAMETER',7X,'PARAMETER',9X,'GRADIENT',
     1       10X,'LAGRANGE',8X,A8/
     2       14X,'NUMBER',10X,'VALUE',28X,'MULTIPLIER',7X,A9/)
 1030 FORMAT(//' ITERATION'//11X,'PHI = VALUE OF OBJECTIVE FUNCTION',
     1        10X,'GNO = NORM OF GRADIENT')
 1040 FORMAT (15X,I3,4X,1PD17.5,1PD17.5,1PD17.5,1PD17.5)
 1050 FORMAT (15X,I3,4X,1PD17.5,1PD17.5,1PD17.5)

      END

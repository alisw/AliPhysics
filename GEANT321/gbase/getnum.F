*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GETNUM(LIST,N)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to count the number of non blank elements        *
C.    *      in the array LIST before the first blank one              *
C.    *                                                                *
C.    *    ==>Called by : GFFGO                                        *
C.    *    ==>Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION LIST(*)
      SAVE IFIRST,IDFLT
      DATA IFIRST/0/
C.    ------------------------------------------------------------------
C.
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH('    ',IDFLT,4,4)
      ENDIF
C
      N = 0
      DO 10 I=1,20
         IF(LIST(I).EQ.IDFLT) GO TO 99
         N=N+1
  10  CONTINUE
C
  99  RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCOL(ICOL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Set line color to ICOL                                   *
C.    *                                                                *
C.    *       if ICOL > 0 set color temporarily                        *
C.    *       if ICOL < 0 set color permanently                        *
C.    *       if ICOL = 0 restore permanent color value                *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDHEAD, GDRAW, GDSPEC,GDXYZ *
C.    *       Author : P.Zanarini ; S.Giani  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      KCOL=MOD(IABS(ICOL),16)
      IF (ICOL.EQ.0) THEN
         CALL MVBITS(LINATP,16,4,LINATT,16)
         CALL MVBITS(LINATP,16,4,ITXATT,3)
      ELSE IF (ICOL.GT.0) THEN
         CALL MVBITS(KCOL,0,4,LINATT,16)
         CALL MVBITS(KCOL,0,4,ITXATT,3)
      ELSE
         CALL MVBITS(KCOL,0,4,LINATT,16)
         CALL MVBITS(KCOL,0,4,LINATP,16)
         CALL MVBITS(KCOL,0,4,ITXATT,3)
      ENDIF
      END

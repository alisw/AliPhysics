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
      SUBROUTINE GDCOL1(ICOL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Change the line color to ICOL                            *
C.    *                                                                *
C.    *       The color table is :                                     *
C.    *                                                                *
C.    *       1 = White (Black is suggested for paper outputs)         *
C.    *       2 = Red                                                  *
C.    *       3 = Green                                                *
C.    *       4 = Blue                                                 *
C.    *                                                                *
C.    *       5 = Yellow                                               *
C.    *       6 = Cyan                                                 *
C.    *       7 = Magenta                                              *
C.    *       8 = Orange                                               *
C.    *                                                                *
C.    *    ==>Called by : GDRAWP, GDRAWT                               *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
C
      SAVE IOLD
      DATA IOLD/-1/
C.
C.    ------------------------------------------------------------------
C.
      LCOL=MOD(IABS(ICOL),16)
      KCOL=LOOKTB(LCOL)
C
      IF (KCOL.NE.IOLD) THEN
         IOLD=KCOL
         CALL ISPLCI(KCOL)
         CALL ISTXCI(KCOL)
         CALL ISPMCI(KCOL)
         CALL ISFACI(KCOL)
      ENDIF
      END

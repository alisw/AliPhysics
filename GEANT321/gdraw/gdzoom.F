*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDZOOM(ZOOMFU,ZOOMFV,UUZ0,VVZ0,UU0,VV0)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Zoom by the factors ZOOMFU,ZOOMFV                        *
C.    *       around screen coordinates UUZ0,VVZ0                      *
C.    *       and places the resulting picture at UU0,VV0              *
C.    *                                                                *
C.    *       After GDZOOM has been called, any drawing operation      *
C.    *       will actually draw accordingly to last GDZOOM parameters *
C.    *                                                                *
C.    *       If ZOOMFU or ZOOMFV is 0 (or negative) the zoom          *
C.    *       parameters (GZUA,...,GZVC in common GCDRAW) are reset    *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF (ZOOMFU.LE.0.OR.ZOOMFV.LE.0) THEN
         GZUA=1
         GZVA=1
         GZUB=0
         GZVB=0
         GZUC=0
         GZVC=0
      ELSE
         GZUA=GZUA*ZOOMFU
         GZVA=GZVA*ZOOMFV
         GZUB=(GZUB+GZUC)*ZOOMFU
         GZVB=(GZVB+GZVC)*ZOOMFV
         GZUC=UU0-UUZ0*ZOOMFU
         GZVC=VV0-VVZ0*ZOOMFV
      ENDIF
C
      END

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
      SUBROUTINE GDXZOO(ISEL1,ZOOMFU,ZOOMFV,UUZ0,VVZ0,UU0,VV0)
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
C.    *       It is now possible to show in different windows the      *
C.    *       `zoomed-views' obtained just clicking in different       *
C.    *       points of the main picture or of the zoomed-pictures     *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author: S. Giani                                         *
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
        II=0
 15     CONTINUE
        II=II+1
        GZUA=1
        GZVA=1
        GZUB=0
        GZVB=0
        GZUC=0
        GZVC=0
        CALL IRQLC(1,2,ISTAT,NT,PX,PY)
        IF(ISTAT.EQ.0)THEN
         DO 16 J=2,II
          CALL ICLWK(J)
 16      CONTINUE
         RETURN
        ENDIF
        IF(ISTAT.EQ.1)THEN
         CALL IDAWK(1)
 17      CONTINUE
         UUZ0=PX
         VVZ0=PY
         GZUA=GZUA*ZOOMFU
         GZVA=GZVA*ZOOMFV
         GZUB=(GZUB+GZUC)*ZOOMFU
         GZVB=(GZVB+GZVC)*ZOOMFV
         GZUC=UU0-UUZ0*ZOOMFU
         GZVC=VV0-VVZ0*ZOOMFV
         CALL IOPWK(II+1,1,II+1)
         CALL IACWK(II+1)
         CALL GDSHOW(ISEL1)
         CALL IRQLC(II+1,2,ISTA2,NT,PX,PY)
         IF(ISTA2.EQ.0)THEN
          CALL IDAWK(II+1)
          CALL IACWK(1)
          GOTO 15
         ELSE
          CALL IDAWK(II+1)
          II=II+1
          GOTO 17
         ENDIF
        ENDIF
      ENDIF
C
      END

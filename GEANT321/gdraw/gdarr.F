*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GDARR(X1,Y1,X2,Y2,SIZE,ITEXT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a segment from (X1,Y1) to (X2,Y2)                   *
C.    *       and an arrow (of size SIZE) at the ending point (X2,Y2)  *
C.    *       together with the one-character text ITEXT               *
C.    *                                                                *
C.    *    ==>Called by : GDAXIS                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION U(3),V(3)
      CHARACTER*(*) ITEXT
      SAVE TGAR
      DATA TGAR/0.6/
C.
C.    ------------------------------------------------------------------
C.
      U(1)=X1
      U(2)=X2
      V(1)=Y1
      V(2)=Y2
      XL=SQRT((X2-X1)**2+(Y2-Y1)**2)
C
C             If axis size is too small, then return without drawing it
C
      IF (XL.LE.0.1) GO TO 999
C
C             Draw the segment
C
      CALL GDRAWV(U,V,2)
C
      FSIN=(Y2-Y1)/XL
      FCOS=(X2-X1)/XL
      U(1)=X1+(XL-SIZE)*FCOS-SIZE*TGAR*FSIN
      U(3)=X1+(XL-SIZE)*FCOS+SIZE*TGAR*FSIN
      V(1)=Y1+(XL-SIZE)*FSIN+SIZE*TGAR*FCOS
      V(3)=Y1+(XL-SIZE)*FSIN-SIZE*TGAR*FCOS
      U(2)=X2
      V(2)=Y2
C
C             Draw the arrow
C
      CALL GDRAWV(U,V,3)
C
      X0=X1+(XL+2*SIZE)*FCOS
      Y0=Y1+(XL+2*SIZE)*FSIN
      ANGLE=0.0
      LWIDTH=1
      IOPT=0
C
C             Draw the text
C
      CALL GDRAWT(X0,Y0,ITEXT,SIZE,ANGLE,LWIDTH,IOPT)
C
  999 RETURN
      END

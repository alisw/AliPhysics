*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRAWC(NAME,IAX,CUTVAL,U0,V0,SU,SV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine draws the intersection of the object        *
C.    *       called NAME with the cut plane normal to one of the      *
C.    *       main axes (IAX=1 is X, IAX=2 is Y, IAX=3 is Z)           *
C.    *       and placed at the distance CUTVAL from the origin.       *
C.    *                                                                *
C.    *       The resulting picture is seen from the axis choosen      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDSPEC                      *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      CHARACTER*4 NAME
      DIMENSION THE(3),PHI(3)
      SAVE THE,PHI
      DATA THE/90.,90.,0./
      DATA PHI/180.,90.,0./
C.
C.    ------------------------------------------------------------------
C.
C             Allocate working space
C
      CALL GWORK(NWCUT)
C
      ICUT=IAX
      DCUT=CUTVAL
C
      NSURF=0
      ISURF=1
C
C             ICUTFL=1 informs GDRAW to call GDRWSC instead of GDRAWS
C             (so that we'll have stored surfaces instead of drawn lines)
C             and to use the line buffer logic of GDRAWV
C
C             PSI is fixed to 0.
C
      ICUTFL=1
      CALL GDRAW(NAME,THE(ICUT),PHI(ICUT),0.,U0,V0,SU,SV)
      ICUTFL=0
      END

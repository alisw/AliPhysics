*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:40  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      FUNCTION GARNDM(DUMMY)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       GARNDM = -LOG(RNDM(DUMMY))                               *
C.    *                                                                *
C.    *    ==>Called by : GBREME,GBREMM,GCOMP,GDRAY,GHADRO,GHEISH,     *
C.    *                   GLANDO,GPAIRG,GPAIRM,GRAYL,GLTRAC            *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION RNDM(1)
C.
C.    ------------------------------------------------------------------
C.
      CALL GRNDM(RNDM,1)
      GARNDM=-LOG(RNDM(1))
      END

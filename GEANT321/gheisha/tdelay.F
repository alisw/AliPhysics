*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:59  cernlib
* Geant
*
*
*        SAVED BY F13FES  AT 01/02/84   12:19
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE TDELAY(X)
C
C *** TIME DELAY FOR NUCLEAR REACTIONS ***
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (01-FEB-1984)
C
#include "geant321/mxgkgh.inc"
#include "geant321/s_curpar.inc"
#include "geant321/s_result.inc"
C
      X=0.
      IF(ATNO2.LT.1.5) RETURN
      IF(ATNO2.GT.230.) RETURN
      IF(EK.GT.0.2) RETURN
      X=500.*EXP(-EK/0.04)
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE SELFAB(SPROB)
C
C *** SELF-ABSORBTION IN HEAVY MOLECULES ***
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (11-OCT-1987)
C
#include "geant321/mxgkgh.inc"
#include "geant321/s_curpar.inc"
#include "geant321/s_blank.inc"
C
      SPROB=0.
      EKW=PV(4,MXGKPV)-ABS(PV(5,MXGKPV))
      IF(EKW.LT.5.) RETURN
      ALEKW=LOG(EKW-4.)
      SPROB=0.6*ALEKW
      IF(SPROB.GT.1.) SPROB=1.
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      SUBROUTINE HIGHAB(SPROB)
C
C *** SELF-ABSORBTION IN HEAVY MOLECULES ***
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (11-OCT-1987)
C
C
#include "geant321/mxgkgh.inc"
#include "geant321/s_curpar.inc"
#include "geant321/s_blank.inc"
C
      SPROB=0.
      EKW=PV(4,MXGKPV)-ABS(PV(5,MXGKPV))
      IF(EKW.LT.5.) RETURN
      ALEKW=LOG(EKW-4.)
      BLEKW=0.25-0.02*LOG(ATNO2)
      SPROB=BLEKW*ALEKW
      IF(SPROB.GT.1.) SPROB=1.
      RETURN
      END

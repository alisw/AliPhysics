*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GANNIR
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates positron annihilation in rest                       *
C.    *                                                                *
C.    *    ==>Called by : GTELEC                                       *
C.    *       Author    L.Urban, M.Maire *********                     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
      DIMENSION RNDM(2)
C.
C.    ------------------------------------------------------------------
C.
      KCASE = NAMEC(11)
      IF((IANNI.NE.1).OR.(EMASS.LE.CUTGAM)) THEN
          ISTOP = 2
          DESTEP = DESTEP + 2*EMASS
          RETURN
      ENDIF
         ISTOP = 1
         CALL GRNDM(RNDM,2)
         CT=-1.+2.*RNDM(1)
         ST=SQRT(ABS(1.-CT*CT))
         PHI=TWOPI*RNDM(2)
         NGKINE=NGKINE+1
         GKIN(1,NGKINE)=EMASS*ST*COS(PHI)
         GKIN(2,NGKINE)=EMASS*ST*SIN(PHI)
         GKIN(3,NGKINE)=EMASS*CT
         GKIN(4,NGKINE)=EMASS
         GKIN(5,NGKINE)=1.
         TOFD(NGKINE)=0.
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
         NGKINE=NGKINE+1
         DO 30 J=1,3
  30     GKIN(J,NGKINE)=-GKIN(J,NGKINE-1)
         GKIN(4,NGKINE)=EMASS
         TOFD(NGKINE)=0.
         GKIN(5,NGKINE)=1.
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
      END

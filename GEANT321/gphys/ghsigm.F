*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:36  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      FUNCTION GHSIGM(P,IPART,A)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Returns absorption cross-section in millibarns               *
C.    *   for a particle with                                          *
C.    *                momentum P (GEV/C)                              *
C.    *                GEANT type IPART                                *
C.    *      on a nucleus of atomic number A                           *
C.    *                                                                *
C.    *    ==>Called by : GPHADR                                       *
C.    *       Author    R.Barlow  *********                            *
C.    *       Modified to deal with K0 and extend to 41 particle       *
C.    *                types by John Allison, 31/1/84. *               *
C.    *                                                                *
C.    ******************************************************************
C.
C              The array NTYPE maps the GEANT particle type (IPART)
C              onto GHSIG's internal particle type.  Simple assumptions
C              are made for highly flavoured particles.  The K0 is
C              given a non-existent internal type 6 which causes the
C              mean of the K+ and K- cross sections to be calculated.
C
      DIMENSION NTYPE(41)
      DATA NTYPE/6*0,3*5,6,3,4,2*1,2,6,
     +  8*1,8*2,2*0,4,3,2*6,2*3,1/
C.
C.    ------------------------------------------------------------------
C.
      GHSIGM = 1.E-20
      IF (IPART .GT. 41) RETURN
      ITYPE  = NTYPE(IPART)
      IF (ITYPE .EQ. 0) RETURN
C
      IF (ITYPE .NE. 6)THEN
         GHSIGM = GHSIG(P, ITYPE, A)
      ELSE
         GHSIGM = 0.5 * (GHSIG(P, 3, A) + GHSIG(P, 4, A))
      ENDIF
C
      END

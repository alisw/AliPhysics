*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      FUNCTION GPHSG2(E)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       FUNCTION TO COMPUTE THE PHOTOEFFECT TOTAL CROSS-SECTION  *
C.    *       AS A FUNCTION OF Z AND E                                 *
C.    *                                                                *
C.    *    ==>CALLED BY : special dE/dx                                *
C.    *                                                                *
C.    *    (1) Sandia parametrizsation is used;                        *
C.    *    (2) GPHSG2 in Megabarns                                     *
C.    *    (3) E in GeV                                                *
C.    *                                                                *
C.    *     J. Chwastowski 31.03.93                                    *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#include "geant321/gc10ev.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION EINV,ECUR,ONE,RES,ZERO
#endif
      PARAMETER (ONE=1,ZERO=0)
C.
C.    ------------------------------------------------------------------
C.
      RES = 0.
      IF(E.LT.G10EV) GO TO 30
C Use Sandia data
      JPHXS  = LQ(JPHOT-1)
      NZ     = Q(JPHXS+1)
      ECUR   = E*1E6
      EINV   = ONE/ECUR
      JWEIGH = JPHXS+2*NZ+1
      DO 20 JZ=1,NZ
         JPHFN   = LQ(JPHXS-JZ)
         IPOINT  = JPHFN+1
         IMAX    = Q(IPOINT)
         IPOINT  = IPOINT-4
         DO 10 I = 1,IMAX
            IPOINT = IPOINT+5
            IF(ECUR.LT.Q(IPOINT)) THEN
               J   = IPOINT+1
               RES = RES+Q(JWEIGH+JZ)*
     +         EINV*(Q(J)+EINV*(Q(J+1)+EINV*(Q(J+2)+EINV*Q(J+3))))
               GO TO 20
            ENDIF
   10    CONTINUE
   20 CONTINUE
C RES value is in cm**2/g
C Now in Megabarns
   30 GPHSG2 = MAX(RES,ZERO)*A/AVO*1E-6
C
      END

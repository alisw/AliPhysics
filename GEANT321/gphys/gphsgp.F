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
      FUNCTION GPHSGP(IELEM,E)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       FUNCTION TO COMPUTE THE PHOTOEFFECT TOTAL CROSS-SECTION  *
C.    *       AS A FUNCTION OF Z AND E                                 *
C.    *                                                                *
C.    *    ==>CALLED BY : special dE/dx                                *
C.    *                                                                *
C.    *    (1) Sandia parametrizsation is used;                        *
C.    *    (2) GPHSGP in cm**/g                                        *
C.    *    (3) E in GeV                                                *
C.    *                                                                *
C.    *     J. Chwastowski 31.03.93                                    *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gc10ev.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION EINV,ECUR,ONE
#endif
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      RES = 0.
      IF(E.LT.G10EV) GO TO 20
C Use Sandia data
      JPHXS  = LQ(JPHOT-1)
      IF(IELEM.GT.IQ(JPHXS-2)) THEN
         WRITE(CHMAIL,10010) IELEM, IQ(JPHXS-2)
10010 FORMAT(' *** GPHSGP : element N ',I6,' requested, ',
     +'but only ',I6,' elements in present mixture')
         CALL GMAIL(0,0)
         GOTO 20
      ENDIF
      JPHXSI = LQ(JPHXS-IELEM)
      IPOINT = JPHXSI+1
      IMAX = Q(IPOINT)
      IPOINT = IPOINT+1
      ECUR = E*1.E6
      IF(ECUR.LT.Q(IPOINT)) GO TO 20
      EINV = ONE/ECUR
      DO 10 I = 2,IMAX
         IPOINT = IPOINT+5
         IF(ECUR.LT.Q(IPOINT)) THEN
            J = IPOINT+1
            RES = EINV*(Q(J)+EINV*(Q(J+1)+EINV*(Q(J+2)+EINV*Q(J+3))))
            GO TO 20
         ENDIF
   10 CONTINUE
C This value is in cm**2/g
   20 GPHSGP = MAX(RES,0.)
C
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCURS(U0,V0,MCHAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       When this routine is called in interactive               *
C.    *       version of GEANT3 and with PIGS or GKS option set,       *
C.    *       it draws on screen the present position of the cursor,   *
C.    *       then moving UP-DOWN and LEFT-RIGHT arrows (or similars)  *
C.    *       the user goes on the screen point he wants be selected   *
C.    *       and press a key (but not <CR>) to hold that point.       *
C.    *       If the user types <CR>, it has to press a key again.     *
C.    *       The routine returns then in U0,V0 the user coordinates   *
C.    *       of the graphics cursor, and in MCHAR the ascii value     *
C.    *       of the key pressed.                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDTREE                      *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C
C.    ------------------------------------------------------------------
C.
*
* Request the locator position
*
      CALL IRQLC(1,1,ISTAT,ITNR,U0,V0)
      MCHAR=ISTAT
      END

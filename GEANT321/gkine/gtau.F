*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GTAU
C.
C.    ******************************************************************
C.    *                                                                *
C.    *            Generates the process e+e- -> tau+ tau-             *
C.    *             can be called from GUKINE                          *
C.    *                                                                *
C.    *       Author    R.Hemingway  26 Sept 1985                      *
C.    *       Obsoleted in 3.15                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
      WRITE(CHMAIL,10000)
10000 FORMAT('Obsolete GTAU reached. This routine is now dummy.')
      CALL GMAIL(1,1)
C
      END

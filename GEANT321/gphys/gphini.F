*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GPHINI
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Initialize material constants for                        *
C.    *       the photoelectric effect                                 *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author    J. Chwastowski                                 *
C.    *                                                                *
C.    ******************************************************************
C.
      CALL GPHRIN
      CALL GSHLIN
      CALL GPHXIN
      END

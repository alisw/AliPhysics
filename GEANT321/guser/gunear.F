*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUNEAR(ISEARC, ICALL, X, JNEAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    User search                                                 *
C.    *       ISEARC to identify the given volume                      *
C.    *       ICALL  to identify the calling routine                   *
C.    *              1 GMEDIA like                                     *
C.    *              2 GNEXT like                                      *
C.    *       X      coordinates (+direction for ICALL=2)              *
C.    *       JNEAR  address of default list of neighbours             *
C.    *              (list to be overwriten by user)                   *
C.    *                                                                *
C.    *    Called by : GFTRAC, GINVOL, GTMEDI, GTNEXT, GNEXT, GMEDIA   *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(*)
C.    ------------------------------------------------------------------
*                                                             END GUNEAR
      END

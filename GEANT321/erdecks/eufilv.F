*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE EUFILV (N, EIN, CNAMV, NUMV, IOVL)
*
C.    ******************************************************************
C.    *                                                                *
C.    *    User routine to fill the input values of the commons :      *
C.    *               /EROPTS/, /EROPTC/ and /ERTRIO/ for CHOPT = 'V'  *
C.    *         N     Number of predictions where to store results     *
C.    *         EIN   Input error matrix                               *
C.    *        CNAMV  Volume name of the prediction                    *
C.    *        NUMV   Volume number (if 0 = all volumes)               *
C.    *        IOVL   = 1  prediction when entering in the volume      *
C.    *               = 2  prediction when leaving the volume          *
C.    *                                                                *
C.    *    ==>Called by : USER (before calling ERTRAK)                 *
C.    *       Author    M.Maire, E.Nagy  *********                     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/ertrio.inc"
      DIMENSION EIN(15), NUMV(N), IOVL(N)
      CHARACTER * 4 CNAMV(N)
*
      NEPRED = MIN0 (N, MXPRED)
*
      CALL UCOPY (EIN(1), ERRIN(1),  15)
      DO I=1,NEPRED
        CALL UCTOH (CNAMV(I),NAMEER(I),4,4)
      ENDDO
      CALL UCOPY (NUMV(1),NUMVER(1),NEPRED)
      CALL UCOPY (IOVL(1),IOVLER(1),NEPRED)
*
*                                                           END EUFILV
  99  END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:34  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERBCTR (X)
*
C.    ******************************************************************
C.    *                                                                *
C.    *       After the back-tracking transforms                       *
C.    *                  the transport matrix back to the original     *
C.    *                  particle direction                            *
C.    *                                                                *
C.    *       Arguments                                                *
C.    *           X        - Transport Matrix                          *
C.    *                                                                *
C.    *    ==>Called by : ERPEND                                       *
C.    *       Author    E.Nagy  *********                              *
C.    *                                                                *
C.    ******************************************************************
*
         DOUBLE PRECISION X(5,5)
*
      DO 10 I = 1, 5
         DO 5 K = 1, 5
            IF (I.EQ.2) X(I,K) = -X(I,K)
            IF (K.EQ.2) X(I,K) = -X(I,K)
            IF (I.EQ.4) X(I,K) = -X(I,K)
            IF (K.EQ.4) X(I,K) = -X(I,K)
    5    CONTINUE
   10 CONTINUE
*                                                           END ERBCTR
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGTSET(N,T,IREP)
************************************************************************
*                                                                      *
*     Name: CGTSET                                                     *
*     Author: E. Chernyaev                       Date:    31.07.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Set transformation to screen coordinates and           *
*               find back transformation to world coordinates          *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: N - transformation number                                 *
*            T(4,3) - transformation matrix                            *
*                                                                      *
*     Output: IREP - reply ( 0 - O.K.                                  *
*                           -1 - if illigal N,                         *
*                           -2 - if DET(T)=0.)                         *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgctra.inc"
      REAL T(4,3)
*-          P A R A M E T E R S   C H E C K
      IREP   = 0
      IF (N.LE.0 .OR. N.GT.NTMAX)       IREP = -1
      DET = T(1,1)*T(2,2)*T(3,3) - T(3,1)*T(2,2)*T(1,3) +
     +      T(3,1)*T(1,2)*T(2,3) - T(1,1)*T(3,2)*T(2,3) +
     +      T(2,1)*T(3,2)*T(1,3) - T(2,1)*T(1,2)*T(3,3)
*      IF (ABS(DET) .LT. 1./DELWOR)      IREP = -2
      IF (IREP .LT. 0)                  GOTO 999
*           S E T   T R A N S F O R M A T I O N   M A T R I X
      DO 200 I=1,4
        DO 100 J=1,3
          TSCRN(I,J,N) = T(I,J)
  100     CONTINUE
  200   CONTINUE
*           F I N D   M A T R I X   F O R   B A C K   TRANSFORMATION
      TBACK(1,1,N) = (T(2,2)*T(3,3) - T(2,3)*T(3,2)) / DET
      TBACK(2,1,N) =-(T(2,1)*T(3,3) - T(2,3)*T(3,1)) / DET
      TBACK(3,1,N) = (T(2,1)*T(3,2) - T(2,2)*T(3,1)) / DET
      TBACK(4,1,N) =-(T(2,1)*T(3,2)*T(4,3) - T(4,1)*T(3,2)*T(2,3) +
     +                T(4,1)*T(2,2)*T(3,3) - T(2,1)*T(4,2)*T(3,3) +
     +                T(3,1)*T(4,2)*T(2,3) - T(3,1)*T(2,2)*T(4,3))/DET
      TBACK(1,2,N) =-(T(1,2)*T(3,3) - T(1,3)*T(3,2)) / DET
      TBACK(2,2,N) = (T(1,1)*T(3,3) - T(1,3)*T(3,1)) / DET
      TBACK(3,2,N) =-(T(1,1)*T(3,2) - T(1,2)*T(3,1)) / DET
      TBACK(4,2,N) = (T(1,1)*T(3,2)*T(4,3) - T(4,1)*T(3,2)*T(1,3) +
     +                T(4,1)*T(1,2)*T(3,3) - T(1,1)*T(4,2)*T(3,3) +
     +                T(3,1)*T(4,2)*T(1,3) - T(3,1)*T(1,2)*T(4,3))/DET
      TBACK(1,3,N) = (T(1,2)*T(2,3) - T(1,3)*T(2,2)) / DET
      TBACK(2,3,N) =-(T(1,1)*T(2,3) - T(1,3)*T(2,1)) / DET
      TBACK(3,3,N) = (T(1,1)*T(2,2) - T(1,2)*T(2,1)) / DET
      TBACK(4,3,N) =-(T(1,1)*T(2,2)*T(4,3) - T(4,1)*T(2,2)*T(1,3) +
     +                T(4,1)*T(1,2)*T(2,3) - T(1,1)*T(4,2)*T(2,3) +
     +                T(2,1)*T(4,2)*T(1,3) - T(2,1)*T(1,2)*T(4,3))/DET
*
  999 RETURN
      END

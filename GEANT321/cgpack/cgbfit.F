*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGBFIT(FACE,ABCD,NT)
************************************************************************
*                                                                      *
*     Name: CGBFIT                                                     *
*     Author: E. Chernyaev                       Date:    15.03.89     *
*                                                Revised: S.Giani 1992 *
*                                                                      *
*     Function: Find all intersection points of face with plane        *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: FACE(*) - face                                            *
*            ABCD(4) - plane                                           *
*                                                                      *
*     Output: NT - number of intersections                             *
*                                                                      *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgcedg.inc"
#include "geant321/cgdelt.inc"
      REAL      FACE(*),ABCD(4)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION  A,B,C,D
#endif
      INTEGER   KCG(6)
      DATA      KCG/KCGX1,KCGY1,KCGZ1,KCGX2,KCGY2,KCGZ2/
*-
      ERROR  = EEWOR
      A      = ABCD(1)
      B      = ABCD(2)
      C      = ABCD(3)
      D      = ABCD(4)
*      D1=0
*      D2=0
      XDELT  = FACE(KCGBB)*C - B*FACE(KCGCC)
      YDELT  = FACE(KCGCC)*A - C*FACE(KCGAA)
      ZDELT  = FACE(KCGAA)*B - A*FACE(KCGBB)
      K      = 1
      IF (ABS(YDELT) .GT. ABS(XDELT))           K = 2
      IF (ABS(ZDELT) .GT. ABS(DELTA(K)))        K = 3
      IF(ABS(DELTA(K)).LT.ERROR)RETURN
      DD     = 1./ABS(DELTA(K))
      DELTA(1)  = DELTA(1) * DD
      DELTA(2)  = DELTA(2) * DD
      DELTA(3)  = DELTA(3) * DD
      KCG1   = KCG(K)
      KCG2   = KCG(K+3)
*        F I N D   I N T E R S E C T I O N S   P O I N T S
      J      = LCGFAC - LCGEDG
      NT     = 0
      NEDGE  = FACE(KCGNE)
      DO 600 NE=1,NEDGE
        J      = J + LCGEDG
        D1     = A*FACE(J+KCGX1) + B*FACE(J+KCGY1) + C*FACE(J+KCGZ1) + D
        D2     = A*FACE(J+KCGX2) + B*FACE(J+KCGY2) + C*FACE(J+KCGZ2) + D
        IF (ABS(D1) .LE. ERROR)         D1 = 0.
        IF (ABS(D2) .LE. ERROR)         D2 = 0.
        IF (D1*D2 .GT. 0.)              GOTO 600
        NT     = NT + 1
        IF (D1.EQ.0. .AND. D2.EQ.0.)    GOTO 500
        IF (D1.LT.0. .OR.  D2.GT.0.)    TTT(2,NT) =-1.
        IF (D1.GT.0. .OR.  D2.LT.0.)    TTT(2,NT) =+1.
*           D 1   O R   D 2   I S   N O T   Z E R O
        IF (NT .NE. 1)                  GOTO 400
        XA     = (D2*FACE(J+KCGX1)-D1*FACE(J+KCGX2)) / (D2 - D1)
        YA     = (D2*FACE(J+KCGY1)-D1*FACE(J+KCGY2)) / (D2 - D1)
        ZA     = (D2*FACE(J+KCGZ1)-D1*FACE(J+KCGZ2)) / (D2 - D1)
        TTT(1,1)   = 0.
        GOTO 600
  400   X      = (D2*FACE(J+KCG1)-D1*FACE(J+KCG2)) / (D2 - D1)
        TTT(1,NT)  = (X - AA(K)) / DELTA(K)
        GOTO 600
*           D 1   A N D   D 2   B O T H   A R E   Z E R O
  500   IF (NT .NE. 1)                  GOTO 510
        XA     = FACE(J+KCGX1)
        YA     = FACE(J+KCGY1)
        ZA     = FACE(J+KCGZ1)
        TTT(1,1)   = 0.
        TTT(2,1)   = 0.
        GOTO 520
  510   X        = FACE(J+KCG1)
        TTT(1,NT)  = (X - AA(K)) / DELTA(K)
        TTT(2,NT)  = 0.
  520   NT       = NT + 1
        X        = FACE(J+KCG2)
        TTT(1,NT)  = (X - AA(K)) / DELTA(K)
        TTT(2,NT)  = 2.
        IF (TTT(1,NT) .GT. TTT(1,NT-1)) GOTO 600
        TTT(2,NT-1)= 2.
        TTT(2,NT)  = 0.
  600   CONTINUE
      RETURN
      END

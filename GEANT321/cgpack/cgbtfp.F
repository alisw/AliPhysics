*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGBTFP(FACE,ABCD,CC)
************************************************************************
*                                                                      *
*     Name: CGBTFP                                                     *
*     Author: E. Chernyaev                       Date:    25.10.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Test face against plane                                *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: FACE(*) - face                                            *
*            ABCD(4) - plane                                           *
*                                                                      *
*     Output: C(*) - new edges                                         *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgcedg.inc"
#include "geant321/cgdelt.inc"
      REAL      FACE(*),ABCD(4),CC(*)
*-
      NMAX   = CC(KCGAF)
      NEDGE  = CC(KCGNE)
      JC     = LCGFAC + NEDGE*LCGEDG
      NEDGE  = FACE(KCGNE)
      J      = LCGFAC
*
**          T E S T :   F A C E   L I E   I N   T H E   SAME PLANE
*
      IF (ABS(FACE(KCGAA)-ABCD(1)) .GT. EEWOR)  GOTO 200
      IF (ABS(FACE(KCGBB)-ABCD(2)) .GT. EEWOR)  GOTO 200
      IF (ABS(FACE(KCGCC)-ABCD(3)) .GT. EEWOR)  GOTO 200
      IF (ABS(FACE(KCGDD)-ABCD(4)) .GT. EEWOR)  GOTO 999
      IF (NMAX .LT. JC+NEDGE*LCGEDG)            GOTO 998
      DO 110 NE=1,NEDGE
        CC(JC+KCGAE) =-1.
        IF (ABS(FACE(J+KCGX1)-FACE(J+KCGX2)) .GT. EEWOR) GOTO 101
        IF (ABS(FACE(J+KCGY1)-FACE(J+KCGY2)) .GT. EEWOR) GOTO 101
        IF (ABS(FACE(J+KCGZ1)-FACE(J+KCGZ2)) .GT. EEWOR) GOTO 101
        PRINT *,' CGBTFP: very small face'
  101   CC(JC+KCGX1) = FACE(J+KCGX2)
        CC(JC+KCGY1) = FACE(J+KCGY2)
        CC(JC+KCGZ1) = FACE(J+KCGZ2)
        CC(JC+KCGX2) = FACE(J+KCGX1)
        CC(JC+KCGY2) = FACE(J+KCGY1)
        CC(JC+KCGZ2) = FACE(J+KCGZ1)
        JC     = JC + LCGEDG
        J      = J  + LCGEDG
  110   CONTINUE
      CC(KCGNE) = CC(KCGNE) + NEDGE
      GOTO 999
*
**          T E S T :   F A C E   L I E   I N   T H E   INVERSE PLANE
*
  200 IF (ABS(FACE(KCGAA)+ABCD(1)) .GT. EEWOR)  GOTO 300
      IF (ABS(FACE(KCGBB)+ABCD(2)) .GT. EEWOR)  GOTO 300
      IF (ABS(FACE(KCGCC)+ABCD(3)) .GT. EEWOR)  GOTO 300
      IF (ABS(FACE(KCGDD)+ABCD(4)) .GT. EEWOR)  GOTO 999
      IF (NMAX .LT. JC+NEDGE*LCGEDG)            GOTO 998
      DO 210 NE=1,NEDGE
        IF (ABS(FACE(J+KCGX1)-FACE(J+KCGX2)) .GT. EEWOR) GOTO 201
        IF (ABS(FACE(J+KCGY1)-FACE(J+KCGY2)) .GT. EEWOR) GOTO 201
        IF (ABS(FACE(J+KCGZ1)-FACE(J+KCGZ2)) .GT. EEWOR) GOTO 201
        PRINT *,' CGBTFP: very small face'
  201   CC(JC+KCGAE) = 0.
        CC(JC+KCGX1) = FACE(J+KCGX1)
        CC(JC+KCGY1) = FACE(J+KCGY1)
        CC(JC+KCGZ1) = FACE(J+KCGZ1)
        CC(JC+KCGX2) = FACE(J+KCGX2)
        CC(JC+KCGY2) = FACE(J+KCGY2)
        CC(JC+KCGZ2) = FACE(J+KCGZ2)
        JC     = JC + LCGEDG
        J      = J  + LCGEDG
  210   CONTINUE
      CC(KCGNE) = CC(KCGNE) + NEDGE
      GOTO 999
*           M A K E   F A C E   S L I C I N G
  300 CALL CGBFIT(FACE,ABCD,NT)
      IF (NT .EQ. 0)                    GOTO 999
*           P R E P A R E   E D G E S
      TMIN  = +1.
      TMAX  = -1.
      CALL CGBTTT('LE',TMIN,TMAX,NT,NEDGE)
      IF (NEDGE .EQ. 0)                 GOTO 999
      DO 700 NE=1,NEDGE
        IF (NMAX .LT. JC+LCGEDG)        GOTO 998
        CC(KCGNE)    = CC(KCGNE) + 1.
        CC(JC+KCGAE) = 0.
        CC(JC+KCGX1) = XA + XDELT*TTT(1,NE)
        CC(JC+KCGY1) = YA + YDELT*TTT(1,NE)
        CC(JC+KCGZ1) = ZA + ZDELT*TTT(1,NE)
        CC(JC+KCGX2) = XA + XDELT*TTT(2,NE)
        CC(JC+KCGY2) = YA + YDELT*TTT(2,NE)
        CC(JC+KCGZ2) = ZA + ZDELT*TTT(2,NE)
        JC     = JC + LCGEDG
  700   CONTINUE
      GOTO 999
*
  998 CC(KCGAF) =-1.
  999 RETURN
      END

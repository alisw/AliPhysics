*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGSHEL(ATRINV,ERROR,NEDGE,EDGE)
************************************************************************
*                                                                      *
*     Name: CGSHEL                                                     *
*     Author: E. Chernyaev                       Date:    13.04.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Shell sort of edges                                    *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input:  ATRINV - atribute decrement for inverse edges            *
*             ERROR  - presision error                                 *
*     Output: NEDGE - number of edges                                  *
*             EDGE(*,*)  - edges                                       *
*                                                                      *
*     Output:                                                          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      EDGE(LCGEDG,*),SAVE(LCGEDG)
*-
      IF (NEDGE .LE. 1)                 GOTO 999
      KE     = 0
      DO 200 NE=1,NEDGE
        IF (EDGE(KCGX1,NE) .GT. EDGE(KCGX2,NE)+ERROR)   GOTO 120
        IF (EDGE(KCGX1,NE) .LT. EDGE(KCGX2,NE)-ERROR)   GOTO 110
        IF (EDGE(KCGY1,NE) .GT. EDGE(KCGY2,NE)+ERROR)   GOTO 120
        IF (EDGE(KCGY1,NE) .LT. EDGE(KCGY2,NE)-ERROR)   GOTO 110
        IF (EDGE(KCGZ1,NE) .GT. EDGE(KCGZ2,NE)+ERROR)   GOTO 120
        IF (EDGE(KCGZ1,NE) .LT. EDGE(KCGZ2,NE)-ERROR)   GOTO 110
        GOTO 200
  110   KE     = KE + 1
        IF (KE .EQ. NE)         GOTO 200
        EDGE(KCGAE,KE) = EDGE(KCGAE,NE)
        EDGE(KCGX1,KE) = EDGE(KCGX1,NE)
        EDGE(KCGY1,KE) = EDGE(KCGY1,NE)
        EDGE(KCGZ1,KE) = EDGE(KCGZ1,NE)
        EDGE(KCGX2,KE) = EDGE(KCGX2,NE)
        EDGE(KCGY2,KE) = EDGE(KCGY2,NE)
        EDGE(KCGZ2,KE) = EDGE(KCGZ2,NE)
        GOTO 200
  120   KE     = KE + 1
        EDGE(KCGAE,KE) = EDGE(KCGAE,NE)-ATRINV
        X      = EDGE(KCGX1,NE)
        Y      = EDGE(KCGY1,NE)
        Z      = EDGE(KCGZ1,NE)
        EDGE(KCGX1,KE) = EDGE(KCGX2,NE)
        EDGE(KCGY1,KE) = EDGE(KCGY2,NE)
        EDGE(KCGZ1,KE) = EDGE(KCGZ2,NE)
        EDGE(KCGX2,KE) = X
        EDGE(KCGY2,KE) = Y
        EDGE(KCGZ2,KE) = Z
  200   CONTINUE
      NEDGE = KE
*
**          S H E L L   S O R T   O F   E D G E S
*
      IF (NEDGE .LE. 1)                 GOTO 999
      ISTEP  = 1
  210 ISTEP  = ISTEP*3 + 1
      IF (ISTEP*2 .LT. NEDGE)           GOTO 210
*
  300 ISTEP  = ISTEP/3
      DO 500 I=1,NEDGE-ISTEP
        J1    = I
        J2    = I + ISTEP
*           I F  (E D G E (J 1)  .L E.  E D G E (J 2))  G O T O  5 0 0
        IF (EDGE(KCGX1,J1) .LT. EDGE(KCGX1,J2)-ERROR)   GOTO 500
        IF (EDGE(KCGX1,J1) .GT. EDGE(KCGX1,J2)+ERROR)   GOTO 350
        IF (EDGE(KCGY1,J1) .LT. EDGE(KCGY1,J2)-ERROR)   GOTO 500
        IF (EDGE(KCGY1,J1) .GT. EDGE(KCGY1,J2)+ERROR)   GOTO 350
        IF (EDGE(KCGZ1,J1) .LT. EDGE(KCGZ1,J2)-ERROR)   GOTO 500
        IF (EDGE(KCGZ1,J1) .GT. EDGE(KCGZ1,J2)+ERROR)   GOTO 350
        IF (EDGE(KCGX2,J1) .LT. EDGE(KCGX2,J2)-ERROR)   GOTO 500
        IF (EDGE(KCGX2,J1) .GT. EDGE(KCGX2,J2)+ERROR)   GOTO 350
        IF (EDGE(KCGY2,J1) .LT. EDGE(KCGY2,J2)-ERROR)   GOTO 500
        IF (EDGE(KCGY2,J1) .GT. EDGE(KCGY2,J2)+ERROR)   GOTO 350
        IF (EDGE(KCGZ2,J1) .LT. EDGE(KCGZ2,J2)-ERROR)   GOTO 500
        IF (EDGE(KCGZ2,J1) .GT. EDGE(KCGZ2,J2)+ERROR)   GOTO 350
        GOTO 500
*            S A V E = E D G E (J 2)
  350   SAVE(KCGAE) = EDGE(KCGAE,J2)
        SAVE(KCGX1) = EDGE(KCGX1,J2)
        SAVE(KCGY1) = EDGE(KCGY1,J2)
        SAVE(KCGZ1) = EDGE(KCGZ1,J2)
        SAVE(KCGX2) = EDGE(KCGX2,J2)
        SAVE(KCGY2) = EDGE(KCGY2,J2)
        SAVE(KCGZ2) = EDGE(KCGZ2,J2)
*            E D G E (J 2) = E D G E (J 1)
  400   EDGE(KCGAE,J2) = EDGE(KCGAE,J1)
        EDGE(KCGX1,J2) = EDGE(KCGX1,J1)
        EDGE(KCGY1,J2) = EDGE(KCGY1,J1)
        EDGE(KCGZ1,J2) = EDGE(KCGZ1,J1)
        EDGE(KCGX2,J2) = EDGE(KCGX2,J1)
        EDGE(KCGY2,J2) = EDGE(KCGY2,J1)
        EDGE(KCGZ2,J2) = EDGE(KCGZ2,J1)
        J2     = J1
        J1     = J1 - ISTEP
        IF (J1 .LE. 0)     GOTO 450
*           I F  (E D G E (J 1) .G T. S A V E)  G O T O  4 0 0
        IF (EDGE(KCGX1,J1) .LT. SAVE(KCGX1)-ERROR)      GOTO 450
        IF (EDGE(KCGX1,J1) .GT. SAVE(KCGX1)+ERROR)      GOTO 400
        IF (EDGE(KCGY1,J1) .LT. SAVE(KCGY1)-ERROR)      GOTO 450
        IF (EDGE(KCGY1,J1) .GT. SAVE(KCGY1)+ERROR)      GOTO 400
        IF (EDGE(KCGZ1,J1) .LT. SAVE(KCGZ1)-ERROR)      GOTO 450
        IF (EDGE(KCGZ1,J1) .GT. SAVE(KCGZ1)+ERROR)      GOTO 400
        IF (EDGE(KCGX2,J1) .LT. SAVE(KCGX2)-ERROR)      GOTO 450
        IF (EDGE(KCGX2,J1) .GT. SAVE(KCGX2)+ERROR)      GOTO 400
        IF (EDGE(KCGY2,J1) .LT. SAVE(KCGY2)-ERROR)      GOTO 450
        IF (EDGE(KCGY2,J1) .GT. SAVE(KCGY2)+ERROR)      GOTO 400
        IF (EDGE(KCGZ2,J1) .LT. SAVE(KCGZ2)-ERROR)      GOTO 450
        IF (EDGE(KCGZ2,J1) .GT. SAVE(KCGZ2)+ERROR)      GOTO 400
*           E D G E (J 2) = S A V E
  450   EDGE(KCGAE,J2) = SAVE(KCGAE)
        EDGE(KCGX1,J2) = SAVE(KCGX1)
        EDGE(KCGY1,J2) = SAVE(KCGY1)
        EDGE(KCGZ1,J2) = SAVE(KCGZ1)
        EDGE(KCGX2,J2) = SAVE(KCGX2)
        EDGE(KCGY2,J2) = SAVE(KCGY2)
        EDGE(KCGZ2,J2) = SAVE(KCGZ2)
  500   CONTINUE
      IF (ISTEP .NE. 1)         GOTO 300
*
  999 RETURN
      END

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
      SUBROUTINE CGBSOR(NEDGE,EDGE)
************************************************************************
*                                                                      *
*     Name: CGBSOR                                                     *
*     Author: E. Chernyaev                       Date:    15.03.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Shell sort of edges                                    *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input:  NEDGE - number of edges                                  *
*     Output: EDGE(*,*)  - edges                                       *
*                                                                      *
*     Output:                                                          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL      EDGE(LCGEDG,*)
*-
      ERROR  = 1.5*EEWOR
      ATRINV = 4.
      CALL CGSHEL(ATRINV,ERROR,NEDGE,EDGE)
*
**          D E L E T E   D O U B L E   E D G E S
**          A T R I B U T E   A N A L I S A T I O N
*
      IF (NEDGE .LE. 0)                 GOTO 999
      IF (NEDGE .EQ. 1)                 GOTO 810
      J       = 1
      NSAME = 1
      DO 800 NE=2,NEDGE
        IF (J .EQ. 0)                   GOTO 720
        IF (ABS(EDGE(KCGX1,NE)-EDGE(KCGX1,J)) .GT. ERROR) GOTO 720
        IF (ABS(EDGE(KCGY1,NE)-EDGE(KCGY1,J)) .GT. ERROR) GOTO 720
        IF (ABS(EDGE(KCGZ1,NE)-EDGE(KCGZ1,J)) .GT. ERROR) GOTO 720
        IF (ABS(EDGE(KCGX2,NE)-EDGE(KCGX2,J)) .GT. ERROR) GOTO 720
        IF (ABS(EDGE(KCGY2,NE)-EDGE(KCGY2,J)) .GT. ERROR) GOTO 720
        IF (ABS(EDGE(KCGZ2,NE)-EDGE(KCGZ2,J)) .GT. ERROR) GOTO 720
*           A T R I B U T E   A N A L I S A T I O N
        AJ     = EDGE(KCGAE,J)
        ANE    = EDGE(KCGAE,NE)
        IF (AJ .EQ. ANE)                GOTO 705
        IF (AJ.LE.-4. .AND. ANE.GT.-4.) GOTO 710
        IF (AJ.GT.-4. .AND. ANE.LE.-4.) GOTO 710
        IF (AJ .GT. ANE)                EDGE(KCGAE,J) = EDGE(KCGAE,NE)
  705   EDGE(KCGX1,J)=((EDGE(KCGX1,J)*NSAME)+EDGE(KCGX1,NE))/(NSAME+1)
        EDGE(KCGY1,J)=((EDGE(KCGY1,J)*NSAME)+EDGE(KCGY1,NE))/(NSAME+1)
        EDGE(KCGZ1,J)=((EDGE(KCGZ1,J)*NSAME)+EDGE(KCGZ1,NE))/(NSAME+1)
        EDGE(KCGX2,J)=((EDGE(KCGX2,J)*NSAME)+EDGE(KCGX2,NE))/(NSAME+1)
        EDGE(KCGY2,J)=((EDGE(KCGY2,J)*NSAME)+EDGE(KCGY2,NE))/(NSAME+1)
        EDGE(KCGZ2,J)=((EDGE(KCGZ2,J)*NSAME)+EDGE(KCGZ2,NE))/(NSAME+1)
        NSAME = NSAME + 1
        GOTO 800
*
  710   J      = J - 1
        GOTO 800
*
  720   J      = J + 1
        EDGE(KCGAE,J) = EDGE(KCGAE,NE)
        EDGE(KCGX1,J) = EDGE(KCGX1,NE)
        EDGE(KCGY1,J) = EDGE(KCGY1,NE)
        EDGE(KCGZ1,J) = EDGE(KCGZ1,NE)
        EDGE(KCGX2,J) = EDGE(KCGX2,NE)
        EDGE(KCGY2,J) = EDGE(KCGY2,NE)
        EDGE(KCGZ2,J) = EDGE(KCGZ2,NE)
        NSAME = 1
  800   CONTINUE
      NEDGE  = J
*
  810 DO 900 NE=1,NEDGE
        ANE    = EDGE(KCGAE,NE)
        IF (ANE .GE. -3.)       GOTO 850
        ANE    = ANE + 4.
        X      = EDGE(KCGX1,NE)
        Y      = EDGE(KCGY1,NE)
        Z      = EDGE(KCGZ1,NE)
        EDGE(KCGX1,NE) = EDGE(KCGX2,NE)
        EDGE(KCGY1,NE) = EDGE(KCGY2,NE)
        EDGE(KCGZ1,NE) = EDGE(KCGZ2,NE)
        EDGE(KCGX2,NE) = X
        EDGE(KCGY2,NE) = Y
        EDGE(KCGZ2,NE) = Z
  850   IF (ANE .LT. -1.)       ANE = ANE + 2.
        EDGE(KCGAE,NE) = ANE
  900   CONTINUE
*
  999 RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE LFITW(X,Y,W,L,KEY,A,B,E)
C
C     TO PERFORM A WEIGHTED STRAIGHT LINE FIT
C
C     FOR FORMULAE USED SEE MENZEL, FORMULAS OF PHYSICS P.116
C
C     FIT IS OF Y=AX+B , WITH S**2 ESTIMATOR E. WEIGHTS ARE IN W.
C     IF KEY=0, POINTS WITH Y=0 ARE IGNORED
C     L IS NO. OF POINTS
C
      DIMENSION X(1),Y(1),W(1)
C
C     CALCULTE SUMS
      IF(L.LE.1) GO TO 1
      W2=0.
      W2X=0.
      W2Y=0.
      W2XY=0.
      W2X2=0.
      W2Y2=0.
      ICNT=0
      DO 2 J=1,L
      IF(Y(J).EQ.0..AND.KEY.EQ.0) GO TO 2
      WW=W(J)*W(J)
      W2=WW+W2
      WWF=WW*X(J)
      W2X=WWF+W2X
      W2X2=WWF*X(J)+W2X2
      W2XY=WWF*Y(J)+W2XY
      WWFI=WW*Y(J)
      W2Y=WWFI+W2Y
      W2Y2=WWFI*Y(J)+W2Y2
      ICNT=ICNT+1
    2 CONTINUE
C
C     FIT PARAMETERS
      A=(W2XY-W2X*W2Y/W2)/(W2X2-W2X**2/W2)
      B=(W2Y-A*W2X)/W2
      IF(ICNT.LE.2) GO TO 3
      E=(W2Y2-W2Y**2/W2-(W2XY-W2X*W2Y/W2)**2/(W2X2-W2X**2/W2))/(ICNT-2)
      GO TO 4
C
C     ISUFFICIENT POINTS
    1 A=0.
      B=0.
    3 E=0.
    4 RETURN
      END

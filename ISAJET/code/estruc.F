#include "isajet/pilot.h"
      FUNCTION ESTRUC(X,QS)
C
C     THIS IS ELECTRON PARTON DISTRIBUTION FUNCTION;
C     SAME AS USED IN PYTHIA; NOTE! ESTRUC=0 FOR X>.999999
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL ESTRUC,AL,PI,AME,QS,X,BT,XM,T,A,B
C
      AL=1./128.
      PI=4*ATAN(1.)
      AME=.511E-3
      BT=2*AL/PI*(LOG(QS/AME/AME)-1.)
C     KLEISS/SJOSTRAND PRESCRIPTION
C      IF (X.LE..9999) THEN
C        ESTRUC=BT/2.*(1.-X)**(BT/2.-1.)
C      ELSE IF (X.LE..999999.AND.X.GT..9999) THEN
C        ESTRUC=100.**(BT/2.)/(100.**(BT/2.)-1.)*BT/2.*
C    $          (1.-X)**(BT/2.-1.)
C      ELSE
C        ESTRUC=0.
C      END IF
C     FADIN-KURAEV/DREES PRESCRIPTION
      XM=.998
      IF(X.GT.XM) THEN
        T = (1.+.375*BT)*(1.-XM)**(BT/2.)
        A = ((1.0-BT/2.)*T
     &        -.25*BT*(1.5-XM*(1.+XM/2.)))/(1.-XM)
     &        +.25*BT*(1.0+XM)
        A = 2*A/(1.-XM)
        B = .5*BT*T/(1.-XM) - .25*BT*(1.+XM) - A*XM
        ESTRUC = A*X+B
      ELSE
        ESTRUC = .5*BT*((1.-X)**(.5*BT-1.)) * (1.+.375*BT)
     &            -.25*BT*(1.+X)
      ENDIF
      RETURN
      END

#include "isajet/pilot.h"
      SUBROUTINE WWKIN(AM1,AM2)
C          WPAIR KINEMATICS, INCLUDING DOUBLE PRECISION CONVERSION FOR
C          32-BIT MACHINES. CONVENTION IS THAT SINGLE PRECISION MASSES
C          AM1,AM2 ARE EXACT.
#include "isajet/itapes.inc"
#include "isajet/wwpar.inc"
#include "isajet/jetpar.inc"
C          BASIC KINEMATICS FROM TWOKIN
      CALL TWOKIN(0.,0.,AM1,AM2)
C          WPAIR KINEMATICS -- JUST A COPY FOR CDC BUT CONSTRUCTS A
C          CONSISTENT SET OF DOUBLE PRECISION VARIABLES FOR 32-BIT
C          MACHINES.
      P3WW(1)=PT(1)*COS(PHI(1))
      P3WW(2)=PT(1)*SIN(PHI(1))
      P3WW(3)=P(1)*CTH(1)
      P3WW(5)=AM1
      P4WW(1)=-P3WW(1)
      P4WW(2)=-P3WW(2)
      P4WW(3)=P(2)*CTH(2)
      P4WW(5)=AM2
#if defined(CERNLIB_SINGLE)
      P3WW(4)=SQRT(P3WW(1)**2+P3WW(2)**2+P3WW(3)**2+P3WW(5)**2)
      P4WW(4)=SQRT(P4WW(1)**2+P4WW(2)**2+P4WW(3)**2+P4WW(5)**2)
#endif
#if defined(CERNLIB_DOUBLE)
      P3WW(4)=DSQRT(P3WW(1)**2+P3WW(2)**2+P3WW(3)**2+P3WW(5)**2)
      P4WW(4)=DSQRT(P4WW(1)**2+P4WW(2)**2+P4WW(3)**2+P4WW(5)**2)
#endif
      P1WW(1)=0.
      P1WW(2)=0.
      P1WW(4)=.5*(P3WW(4)+P3WW(3)+P4WW(4)+P4WW(3))
      P1WW(3)=P1WW(4)
      P2WW(1)=0.
      P2WW(2)=0.
      P2WW(4)=.5*(P3WW(4)-P3WW(3)+P4WW(4)-P4WW(3))
      P2WW(3)=-P2WW(4)
C          INVARIANTS
      SWW=+2.*(P1WW(4)*P2WW(4)-P1WW(3)*P2WW(3))
      TWW=-2.*(P1WW(4)*P3WW(4)-P1WW(3)*P3WW(3))+P3WW(5)**2
      UWW=-2.*(P2WW(4)*P3WW(4)-P2WW(3)*P3WW(3))+P3WW(5)**2
      RETURN
      END

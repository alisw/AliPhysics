*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE SETUP(X1,I1,X2,I2,Y1,J1,Y2,J2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PURPOSE                                                            C
C     TO SET UP THE MAPPING CONSTANTS BETWEEN TWO LINEARLY DEPENDANT   C
C     COORDINATE SYSTEMS. SETUP IS PLANNED TO BE USED IN CONJUNCTION   C
C     WITH THE FUNCTION FREARG.                                        C
C   USAGE                                                              C
C     CALL SETUP(X1,I1,X2,I2,Y1,J1,Y2,J2)                              C
C     THIS STATEMENT SHOULD BE EXECUTED BEFORE FREARG IS CALLED.       C
C   PARAMETERS                                                         C
C     WE HAVE TWO DIFFERENT COORDINATE SYSTEMS, THE (X,Y)-SYSTEM AND   C
C     THE (I,J)-SYSTEM.                                                C
C     X AND I ARE LINEARLY RELATED TO EACH OTHER.                      C
C     Y AND J ARE LINEARLY RELATED TO EACH OTHER.                      C
C     THE MAPPING CONSTANTS ARE DEFINED BY REQUIRING THAT              C
C     X=X1 CORRESPONDS TO I=I1,                                        C
C     X=X2 CORRESPONDS TO I=I2,                                        C
C     Y=Y1 CORRESPONDS TO J=J1,                                        C
C     Y=Y2 CORRESPONDS TO J=J2.                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON /J509C2/CX1,CI1,CX2,CI2,CY1,CJ1,CY2,CJ2
C=======================================================================
      CX1=X1
      CI1=I1
      CX2=X2
      CI2=I2
      CY1=Y1
      CJ1=J1
      CY2=Y2
      CJ2=J2
      RETURN
      END

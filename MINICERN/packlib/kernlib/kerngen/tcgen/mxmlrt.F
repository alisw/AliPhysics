*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE MXMLRT (A,B,C,NI,NJ)
C
C CERN PROGLIB# F110    MXMLRT          .VERSION KERNFOR  2.00  720707
C ORIG. 01/01/64 RKB
C
      DIMENSION A(*),B(*),C(*)
C
C--      ENTRY MXMLRT
C
C--                C = A(I,J) X B(J,J) X A*(J,I)
C--                A* STANDS FOR A-TRANSPOSED
C
C        CALL MXMLRT (A,B,C,NI,NJ)     IS EQUIVALENT TO
C             CALL MXMPY (A,B,X,NI,NJ,NJ)
C             CALL MXMPY1 (X,A,C,NI,NJ,NI)
C
C        OR   CALL MXMPY1 (B,A,Y,NJ,NJ,NI)
C             CALL MXMPY (A,Y,C,NI,NJ,NI)
C
      IPA=1
      JPA=NJ
      GO TO 5
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXMLTR
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXMLTR   (A,B,C,NI,NJ)
#endif
C
C--                C = A*(I,J) X B(J,J) X A(J,I)
C
C        CALL MXMLTR (A,B,C,NI,NJ)     IS EQUIVALENT TO
C             CALL MXMPY2 (A,B,X,NI,NJ,NJ)
C             CALL MXMPY (X,A,C,NI,NJ,NI)
C
C        OR   CALL MXMPY (B,A,Y,NJ,NJ,NI)
C             CALL MXMPY2 (A,Y,C,NI,NJ,NI)
C
      IPA=NI
      JPA=1
C
    5 CONTINUE
      IF (NI.LE.0) RETURN
      IF (NJ.LE.0) RETURN
      IC1=1
      IA1=1
         DO 50 II=1,NI
C
      IC=IC1
         DO 10 KC=1,NI
      C(IC)=0.
   10 IC=IC+1
C
      IB1=1
      JA1=1
         DO 40 JJ=1,NJ
C
      IB=IB1
      IA=IA1
      X=0.
         DO 20 KJ=1,NJ
      X=X + A(IA)*B(IB)
      IA=IA+IPA
      IB=IB+NJ
   20 CONTINUE
C
      JA=JA1
      IC=IC1
         DO 30 KI=1,NI
      C(IC)=C(IC) + X*A(JA)
      IC=IC+1
      JA=JA+JPA
   30 CONTINUE
C
      IB1=IB1+1
      JA1=JA1+IPA
   40 CONTINUE
C
      IC1=IC1+NI
      IA1=IA1+JPA
   50 CONTINUE
      RETURN
      END

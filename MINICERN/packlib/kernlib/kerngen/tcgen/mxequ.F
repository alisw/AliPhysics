*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:52  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if !defined(CERNLIB_TCGEN)
      SUBROUTINE MXEQU(A,B,I,J)
C
C CERN PROGLIB# F109    MXEQU           .VERSION KERNFOR  1.0   680220
C ORIG. 01/01/64 RKB ET AL
C
      DIMENSION A(*),B(*)
C
      IF (I.EQ.0) RETURN
C
C--                TRANSFORMATION OF A-MATRIX
   10 LL=1
         DO 17 L=1,I
      DIAG=-1./A(LL)
      LM=L
C
         DO 12 M=1,I
      IF (L.EQ.M) GO TO 12
      A(LM)=A(LM)*DIAG
   12 LM=LM+I
C
      LPLUS1=L+1
      IF (LPLUS1.GT.I) GO TO 20
      L1L=LL+1
C
         DO 16 L1=LPLUS1,I
      L1M=L1
      LM=L
C
         DO 15 M=1,I
      IF (L.EQ.M) GO TO 14
      A(L1M)=A(L1M)+A(L1L)*A(LM)
   14 L1M=L1M+I
   15 LM=LM+I
C
   16 L1L=L1L+1
C
   17 LL=LL+I+1
C
C
#if defined(CERNLIB_ENTRCDC)
         ENTRY MXEQU1
#endif
#if !defined(CERNLIB_ENTRCDC)
         ENTRY MXEQU1  (A,B,I,J)
#endif
C
C--                TRANSFORMATION OF B-MATRIX
   20 IF (I*J.EQ.0) RETURN
C
      ML=1
         DO 26 L=1,I
C
         DO 26 M=1,J
      MN=M
      LN=L
C
         DO 25 N=1,I
      IF (L.EQ.N) GO TO 24
      B(MN)=B(MN)+B(ML)*A(LN)
   24 MN=MN+J
   25 LN=LN+I
C
   26 ML=ML+1
C
C--                PUTS FINAL RESULT IN B-MATRIX
      NN=1
      MN=1
         DO 36 N=1,I
      DIAG=1./A(NN)
C
         DO 35 M=1,J
      B(MN)=B(MN)*DIAG
   35 MN=MN+1
C
   36 NN=NN+I+1
C
      RETURN
      END
#endif

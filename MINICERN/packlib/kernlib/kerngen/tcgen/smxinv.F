*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if !defined(CERNLIB_TCGEN)
      SUBROUTINE SMXINV (A,NDIM,IFAIL)
C
C CERN PROGLIB# F107    SMXINV          .VERSION KERNFOR  1.0   720503
C ORIG. 03/05/72 CL
C
      DIMENSION A(*),INDEX(100),RI(100)
C
      DATA  TOL / 1.E-12/
C
      IFAIL=0
      N=NDIM
      NP1=N+1
         DO 10 I=1,N
   10 INDEX(I)=1
C
         DO 80 I=1,N
C
C--                FIND PIVOT
      PIVOT=0.
      JJ=1
         DO 20 J=1,N
      IF (INDEX(J).EQ.0) GO TO 19
      ELM=ABS (A(JJ))
      IF (ELM.LE.PIVOT) GO TO 19
      PIVOT=ELM
      K=J
      KK=JJ
   19 JJ=JJ+NP1
   20 CONTINUE
      IF (PIVOT/ABS(A(1)).LT.TOL) GO TO 100
      INDEX(K)=0
      PIVOT=-A(KK)
C
C--                ELIMINATION
      KJ=K
      NP=N
C
         DO 70 J=1,N
      IF (J-K) 34,30,34
C
   30 A(KJ)=1./PIVOT
      RI(J)=0.
      NP=1
      GO TO 70
C
   34 ELM=-A(KJ)
   40 RI(J)=ELM/PIVOT
      IF (ELM.EQ.0.) GO TO 50
C
      JL=J
         DO 45 L=1,J
      A(JL)=A(JL)+ELM*RI(L)
   45 JL=JL+N
C
   50 A(KJ)=RI(J)
C
   70 KJ=KJ+NP
C
   80 CONTINUE
C
C--                CHANGE THE SIGN AND PROVISIONAL FILL-UP
      IJ0=1
      JI0=1
         DO 95 I=1,N
      IJ=IJ0
      JI=JI0
C
         DO 90 J=1,I
      A(IJ)=-A(IJ)
      A(JI)=A(IJ)
      IJ=IJ+N
      JI=JI+1
   90 CONTINUE
C
      IJ0=IJ0+1
      JI0=JI0+N
   95 CONTINUE
      RETURN
C
C--                FAILURE RETURN
  100 IFAIL=1
      RETURN
      END
#endif

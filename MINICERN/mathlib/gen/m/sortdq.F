*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_QMIBMVF)
@PROCESS DIRECTIVE('*VDIR:') VECTOR
#endif
      SUBROUTINE SORTDQ(A,NC,NR,NS)
      DOUBLE PRECISION A(NC,NR),X
#if defined(CERNLIB_CRAY)
      DOUBLE PRECISION TEMP
#endif
      INTEGER LT(20),RT(20),R
C
      NCS = ABS(NS)
      LEVEL=1
      LT(1)=1
      RT(1)=NR
   10 L=LT(LEVEL)
      R=RT(LEVEL)
      LEVEL=LEVEL-1
   20 IF(R.LE.L) IF(LEVEL) 90,90,10
C
C   SUBDIVIDE THE INTERVAL L,R
C     L : LOWER LIMIT OF THE INTERVAL (INPUT)
C     R : UPPER LIMIT OF THE INTERVAL (INPUT)
C     J : UPPER LIMIT OF LOWER SUB-INTERVAL (OUTPUT)
C     I : LOWER LIMIT OF UPPER SUB-INTERVAL (OUTPUT)
C
      I=L
      J=R
      M=(L+R)/2
      X=A(NCS,M)
   30 IF(A(NCS,I).GE.X) GO TO 40
      I=I+1
      GO TO 30
   40 IF(A(NCS,J).LE.X) GO TO 50
      J=J-1
      GO TO 40
C
   50 IF(I.GT.J) GO TO 70
#if defined(CERNLIB_DOUBLE)||defined(CERNLIB_CDC)
      CALL USWOP(A(1,I),A(1,J),2*NC)
#endif
#if defined(CERNLIB_CRAY)
      DO 60 K=1, NC
        TEMP   = A(K,I)
        A(K,I) = A(K,J)
        A(K,J) = TEMP
   60 CONTINUE
#endif
      I=I+1
      J=J-1
      IF(I.LE.J) GO TO 30
C
   70 LEVEL=LEVEL+1
      IF((R-I).GE.(J-L)) GO TO 80
      LT(LEVEL)=L
      RT(LEVEL)=J
      L=I
      GO TO 20
   80 LT(LEVEL)=I
      RT(LEVEL)=R
      R=J
      GO TO 20
   90 IF(NS.GT.0) RETURN
C   REVERSE THE ORDER OF THE ROWS IF NS NEGATIVE
      NRH=NR/2
      DO 110 I=1,NRH
#if (defined(CERNLIB_DOUBLE)||defined(CERNLIB_CDC))&&(!defined(CERNLIB_QMIBMVF))
        CALL USWOP(A(1,I),A(1,NR+1-I),2*NC)
#endif
#if defined(CERNLIB_CRAY)||defined(CERNLIB_QMIBMVF)
        DO 100 K=1, NC
          TEMP = A(K,I)
          A(K,I) = A(K,NR+1-I)
          A(K,NR+1-I) = TEMP
  100   CONTINUE
#endif
  110 CONTINUE
      END

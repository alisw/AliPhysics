*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:50  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMIBMVF)
@PROCESS DIRECTIVE('*VDIR:') VECTOR
#endif
      SUBROUTINE SORTI(A,NC,NR,NS)
C
C CERN PROGLIB# M107    SORTI           .VERSION KERNFOR  4.21  890323
C ORIG. 15/11/88 FCA
C
      INTEGER A(NC,NR), HMIN, HMAX, TEMP
C
      NCS=ABS(NS)
      IF(NCS.EQ.0)  GO TO 999
      IF(NCS.GT.NC) GO TO 999
      IF(NR.LE.1)   GO TO 999
      IF(NS.LE.0)   GO TO 31

C----              Ascending order

      DO 30 J=1,NR-1
#if defined(CERNLIB_QMAPO)
        LMIN = (LVSIMI(A(NCS,J),NR-J+1,NC)-1)/NC+J
#else
        LMIN = J
        HMIN = A(NCS,J)
        DO 10 K=J+1,NR
          IF(HMIN.GT.A(NCS,K)) THEN
            HMIN = A(NCS,K)
            LMIN = K
          ENDIF
   10   CONTINUE
#endif
        IF(LMIN.NE.J) THEN
          DO 25 L=LMIN, J, -1
            IF(A(NCS,L).EQ.A(NCS,J)) THEN
#if defined(CERNLIB_QMAPO)||defined(CERNLIB_QMVAX)||(defined(CERNLIB_QMIBM)&&(!defined(CERNLIB_QMIBMVF)))
              CALL USWOP(A(1,LMIN),A(1,L),NC)
#else
             DO 20 K=1,NC
                TEMP      = A(K,LMIN)
                A(K,LMIN) = A(K,L)
                A(K,L)    = TEMP
   20         CONTINUE
#endif
              LMIN = L
            ENDIF
   25     CONTINUE
        ENDIF
   30 CONTINUE
      GO TO 999

C----              Descending order

   31 DO 60 J=1,NR-1
#if defined(CERNLIB_QMAPO)
        LMAX = (LVSIMX(A(NCS,J),NR-J+1,NC)-1)/NC+J
#else
        LMAX = J
        HMAX = A(NCS,J)
        DO 40 K=J+1,NR
          IF(HMAX.LT.A(NCS,K)) THEN
            HMAX = A(NCS,K)
            LMAX = K
          ENDIF
   40   CONTINUE
#endif
        IF(LMAX.NE.J) THEN
          DO 55 L=LMAX, J, -1
            IF(A(NCS,L).EQ.A(NCS,J)) THEN
#if defined(CERNLIB_QMAPO)||defined(CERNLIB_QMVAX)||(defined(CERNLIB_QMIBM)&&(!defined(CERNLIB_QMIBMVF)))
              CALL USWOP(A(1,LMAX),A(1,L),NC)
#else
              DO 50 K=1,NC
                TEMP      = A(K,LMAX)
                A(K,LMAX) = A(K,L)
                A(K,L)    = TEMP
   50         CONTINUE
#endif
              LMAX = L
            ENDIF
   55     CONTINUE
        ENDIF
   60 CONTINUE
  999 RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE PERMU(IA,N)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='PERMUT')
C
      DIMENSION IA(*)
      PARAMETER (IFD = 12)
      DIMENSION IFCT(0:IFD),IV(IFD+1)
      CHARACTER*80 ERRTXT

      DATA (IFCT(I),I=0,IFD) /1,1,2,6,24,120,720,5040,40320,362880,
     1                        3628800,39916800,479001600/

      IF(N .LE. 0) RETURN
      IF(IA(1) .EQ. 0) THEN
       DO 11 I = 1,N
   11  IA(I)=I
       IF(N .EQ. 1) IA(1)=0
      ELSE
       DO 12 K1 = N,2,-1
       K=K1
       IF(IA(K-1) .LT. IA(K)) GO TO 14
   12  CONTINUE
       IA(1)=0
       RETURN
   14  KN=K+N
       DO 15 L = K,KN/2
       IB=IA(KN-L)
       IA(KN-L)=IA(L)
   15  IA(L)=IB
       DO 16 L1 = K,N
       L=L1
       IF(IA(L) .GT. IA(K-1)) GO TO 17
   16  CONTINUE
   17  IB=IA(K-1)
       IA(K-1)=IA(L)
       IA(L)=IB
      ENDIF
      RETURN

      ENTRY PERMUT(NRP,N,IA)

      IF(N .LE. 0) RETURN
      IF(N .GT. IFD) THEN
       WRITE(ERRTXT,101) N
       CALL MTLPRT(NAME,'V202.1',ERRTXT)
      ELSEIF(NRP .GT. IFCT(N)) THEN
       IA(1)=0
       CALL MTLPRT(NAME,'V202.2',
     +             'PERMUTATION OUTSIDE LEXICON REQUESTED')
      ELSE
       DO 21 I = 1,N
   21  IV(I)=I
       IO=NRP-1
       DO 22 M = N-1,1,-1
       IN=IO/IFCT(M)+1
       IO=MOD(IO,IFCT(M))
       IA(N-M)=IV(IN)
       DO 23 I = IN,M
   23  IV(I)=IV(I+1)
   22  CONTINUE
       IA(N)=IV(1)
      END IF
      RETURN

      ENTRY COMBI(IA,N,J)

      IF(N .LE. 0 .OR. J .LE. 0) RETURN
      IF(J .GT. N) THEN
       WRITE(ERRTXT,103) J,N
       CALL MTLPRT(NAME,'V202.3',ERRTXT)
      ELSEIF(IA(1) .EQ. 0) THEN
       DO 31 I = 1,J
   31  IA(I)=I
       IA(J+1)=0
      ELSE
       DO 32 I1 = 1,N
       I=I1
       IF(IA(I+1) .NE. IA(I)+1) GO TO 33
   32  IA(I)=I
   33  IA(I)=IA(I)+1
       IF(IA(J) .EQ. N+1) IA(1)=0
      ENDIF
      RETURN
  101 FORMAT('N = ',I20,' TOO BIG')
  103 FORMAT('J = ',I5,' > ',I5,' = N IS NOT PERMITTED')
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
*
* c344cod1.inc
*
      X=Z
      Y=-I*Z
      LER=.FALSE.
      IF(X .LE. 0 .AND. Y .EQ. 0) THEN
       WRITE(ERRTXT,101) Z
       CALL MTLPRT(NAME,'C344.1',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. A .AND. A .LT. 1)) THEN
       WRITE(ERRTXT,102) A
       CALL MTLPRT(NAME,'C344.2',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. NMAX .AND. NMAX .LE. 100)) THEN
       WRITE(ERRTXT,103) NMAX
       CALL MTLPRT(NAME,'C344.3',ERRTXT)
       LER=.TRUE.
      END IF
      IF(LER) RETURN
      AZ=ABS(Z)
      YA=ABS(Y)
      ZA=X+I*YA
      EPS=HF*Z10**(-ND)
      DO 5 N = 0,NMAX
      RR(N)=0
    5 CBA(N)=0

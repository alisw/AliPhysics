*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
*
* simpscod.inc
*
      DIMENSION F(0:*)

      IF(N2 .LE. 0 .OR. 2*(N2/2) .NE. N2) THEN
       H=0
       WRITE(ERRTXT,101) N2+1
       CALL MTLPRT(NAME,'D101.1',ERRTXT)
      ELSE
*      S1=0
       S1=F(N2-1)
       S2=0
*      DO 1 N = 1,N2-1,2
*   1  S1=S1+F(N)
*      DO 2 N = 2,N2-2,2
*   2  S2=S2+F(N)

       DO 1 N = 1,N2-3,2
         S1=S1+F(N)
         S2=S2+F(N+1)
    1  CONTINUE
*      S1=S1+F(N2-1)

       S1=S1+S1+S2
       H=(F(0)+F(N2)+S1+S1)*(B-A)/(3*N2)
      ENDIF

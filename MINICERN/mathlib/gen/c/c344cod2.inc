*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
*
* c344cod2.inc
*
      D=TC(8)*ND+TC(9)
      Q=0
      R0=AZ
      IF(NMAX .GT. 0) THEN
       V=0.5*D/NMAX
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 6 J = 2,6
    6   T=V*T+TC(J)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       END IF
       Q=NMAX*T
      END IF
      T=1
      IF(YA .LT. D) THEN
       V=TC(10)*(D-YA)/R0
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 7 J = 2,6
    7   T=V*T+TC(J)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       END IF
      END IF
      NU=1+MAX(Q,TC(11)*R0*T)

      A2=A+A
      MU=-1
    2 MU=MU+1
      AL=1
      DO 3 N = 1,NU
      XN=N
    3 AL=AL*(XN+A2)/(XN+1)
      C=(-I)**NU
      R=0
      S=0
      DO 4 N = NU,1,-1
      XN=N
      A2N=2*(A+XN)
      R=ZA/(A2N-ZA*R)
      IF(N .LE. NMAX) RR(N-1)=R
      AL=AL*(XN+1)/(XN+A2)
      S=R*(A2N*AL*C+S)
    4 C=I*C
      CB(0)=F/(1+S)
      DO 10 N = 0,NMAX-1
   10 CB(N+1)=RR(N)*CB(N)
      IF(Y .LT. 0) THEN
       DO 11 N = 0,NMAX
   11  CB(N)=GCONJG(CB(N))
      END IF
      DO 12 N = 0,NMAX
      IF(ABS(CB(N)-CBA(N)) .GT. EPS*ABS(CB(N))) THEN
       DO 13 M = 0,NMAX
   13  CBA(M)=CB(M)
       NU=NU+5
       IF(MU .LE. 50) GO TO 2
       WRITE(ERRTXT,104) Z,A
       CALL MTLPRT(NAME,'C344.4',ERRTXT)
       RETURN
      END IF
   12 CONTINUE
      RETURN
  101 FORMAT('ILLEGAL ARGUMENT Z = ',1P,2D15.8)
  102 FORMAT('ILLEGAL ORDER A = ',1P,D15.8)
  103 FORMAT('ILLEGAL NMAX = ',I5)
  104 FORMAT('NO CONVERGENCE, Z = ',1P,2D13.6,' A = ',D13.6,
     1       ' TRY SMALLER ND')
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:09  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_IBMVM)
      SUBROUTINE WQBSJA(Z,A,NMAX,ND,CB)

#include "gen/imp128.inc"
#include "gen/defc128.inc"
     +    Z,ZA,CB,I,CBA,RR,F,C,S,R,GCONJG,ZZ
      REAL R0,D,T,Q,U,V,TC(11)
      CHARACTER*80 ERRTXT
      CHARACTER NAME*(*)
      PARAMETER (NAME = 'WQBSJA')
      LOGICAL LER
      DIMENSION CB(0:*),CBA(0:100),RR(0:100)

      PARAMETER (I = (0Q0,1Q0))
      PARAMETER (Z0 = 0, Z1 = 1, HF = Z1/2, Z10 = 10)

      DATA TC / 5.7941 E-5,-1.76148E-3, 2.08645E-2,-1.29013E-1,
     1          8.5777 E-1, 1.0125 E+0, 7.75   E-1, 2.3026 E+0,
     2          1.3863 E+0, 7.3576 E-1, 1.3591 E+0/

      GCONJG(ZZ)=QCONJG(ZZ)
#include "c344cod1.inc"
      F=EXP(YA+I*(A*ATAN2(YA,X)-X))*(HF*AZ)**A/QGAMMA(1+A)
#include "c344cod2.inc"
#endif
*
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE WBSJA(Z,A,NMAX,ND,CB)

#include "gen/imp128.inc"
      REAL R0,D,T,Q,U,V,TC
      LOGICAL LER
      DIMENSION Z(2),CB(2,0:*),CBA(2,0:100),RR(2,0:100)
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'WBSJA')
      DIMENSION TC(11)

      PARAMETER (Z0 = 0, Z1 = 1, Z2 = 2, Z10 = 10, HF = Z1/2)

      DATA TC / 5.7941 E-5,-1.76148E-3, 2.08645E-2,-1.29013E-1,
     1          8.5777 E-1, 1.0125 E+0, 7.75   E-1, 2.3026 E+0,
     2          1.3863 E+0, 7.3576 E-1, 1.3591 E+0/

      X=Z(1)
      Y=Z(2)
      LER=.FALSE.
      IF(X .LE. 0 .AND. Y .EQ. 0) THEN
       WRITE(ERRTXT,101) Z(1),Z(2)
       CALL MTLPRT(NAME,'C344.1',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. A .AND. A .LT. 1)) THEN
       WRITE(ERRTXT,102) A
       CALL MTLPRT(NAME,'C344.2',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. NMAX .AND. NMAX .LT. 100)) THEN
       WRITE(ERRTXT,103) NMAX
       CALL MTLPRT(NAME,'C344.3',ERRTXT)
       LER=.TRUE.
      END IF
      IF(LER) RETURN
      RS=X**2+Y**2
      R1=SQRT(RS)
      YA=ABS(Y)
      DEPS=HF*Z10**(-ND)
      DO 5 N = 0,NMAX
      DO 5 I = 1,2
      RR(I,N)=0
    5 CBA(I,N)=0
      PHI=A*ATAN2(YA,X)-X
      C=EXP(YA)*(HF*R1)**A/DGAMMA(1+A)
      SU1=C*COS(PHI)
      SU2=C*SIN(PHI)
      R0=R1
      D=TC(8)*ND+TC(9)
      Q=0
      IF(NMAX .GT. 0) THEN
       V=D/(2*NMAX)
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 6 I = 2,6
    6   T=V*T+TC(I)
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
        DO 7 I = 2,6
    7   T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       END IF
      END IF
      NU=1+MAX(Q,TC(11)*R0*T)

      H1=A+A
      MU=-1
    2 MU=MU+1
      AL=1
      C1=1
      C2=0
      DO 3 N = 1,NU
      DN=N
      AL=AL*(DN+H1)/(DN+1)
      C=-C1
      C1=C2
    3 C2=C
      R1=0
      R2=0
      S1=0
      S2=0
      DO 4 N = NU,1,-1
      DN=N
      H=2*(A+DN)
      C=1/((H-X*R1+YA*R2)**2+(X*R2+YA*R1)**2)
      R1=(H*X-RS*R1)*C
      R2=(H*YA+RS*R2)*C
      IF(N .LE. NMAX) THEN
       RR(1,N-1)=R1
       RR(2,N-1)=R2
      END IF
      AL=AL*(DN+1)/(DN+H1)
      C=H*AL
      A1=C*C1
      A2=C*C2
      C=C1
      C1=-C2
      C2=C
      AS1=A1+S1
      AS2=A2+S2
      S=R1*AS1-R2*AS2
      S2=R1*AS2+R2*AS1
    4 S1=S
      AS1=1+S1
      C=1/(AS1**2+S2**2)
      CB(1,0)=(SU1*AS1+SU2*S2)*C
      CB(2,0)=(SU2*AS1-SU1*S2)*C
      DO 10 N = 0,NMAX-1
      CB(1,N+1)=RR(1,N)*CB(1,N)-RR(2,N)*CB(2,N)
   10 CB(2,N+1)=RR(1,N)*CB(2,N)+RR(2,N)*CB(1,N)
      IF(Y .LT. 0) THEN
       DO 11 N = 0,NMAX
   11  CB(2,N)=-CB(2,N)
      END IF
      DO 12 N = 0,NMAX
      IF(SQRT((CB(1,N)-CBA(1,N))**2+(CB(2,N)-CBA(2,N))**2)
     1       .GT. DEPS*SQRT(CB(1,N)**2+CB(2,N)**2)) THEN
       DO 13 M = 0,NMAX
       CBA(1,M)=CB(1,M)
   13  CBA(2,M)=CB(2,M)
       NU=NU+5
       IF(MU .LE. 50) GO TO 2
       WRITE(ERRTXT,104) Z(1),Z(2),A
       CALL MTLPRT(NAME,'C344.4',ERRTXT)
       RETURN
      END IF
   12 CONTINUE
      RETURN
  101 FORMAT('ILLEGAL ARGUMENT Z = (',1P,D15.8,', ',D15.8,')')
  102 FORMAT('ILLEGAL ORDER A = ',1P,D15.8)
  103 FORMAT('ILLEGAL NMAX = ',I5)
  104 FORMAT('NO CONVERGENCE, Z = (',1P,D13.6,', ',D13.6,')',
     1       ' A = ',D13.6,' TRY SMALLER ND')
      END
#endif

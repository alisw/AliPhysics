*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:08  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_QUAD)
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE QBSJA(X,A,NMAX,ND,B)
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE DBSJA(X,A,NMAX,ND,B)
#endif

#include "gen/imp128.inc"
      REAL SX,D,T,Q,U,V,TC(11)
      CHARACTER*80 ERRTXT
      CHARACTER NAMEJ*(*),NAMEI*(*)
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEJ = 'DBSJA/QBSJA',
     1           NAMEI = 'DBSIA/QBSIA')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAMEJ = 'DBSJA', NAMEI = 'DBSIA')
#endif
      LOGICAL LJA,LIA,LEV,LER
      DIMENSION B(0:*),BA(0:100),RR(0:100)

      PARAMETER (Z1 = 1, HF = Z1/2, Z10 = 10)

      DATA TC / 5.7941 E-5,-1.76148E-3, 2.08645E-2,-1.29013E-1,
     1          8.5777 E-1, 1.0125 E+0, 7.75   E-1, 2.3026 E+0,
     2          1.3863 E+0, 7.3576 E-1, 1.3591 E+0/

      LJA=.TRUE.
      LIA=.FALSE.
      SGN=-1
      GO TO 9

#if defined(CERNLIB_DOUBLE)
      ENTRY QBSIA(X,A,NMAX,ND,B)
#endif
#if !defined(CERNLIB_DOUBLE)
      ENTRY DBSIA(X,A,NMAX,ND,B)
#endif
      LJA=.FALSE.
      LIA=.TRUE.
      SGN=1

    9 LER=.FALSE.
      IF(X .LE. 0) THEN
       WRITE(ERRTXT,101) X
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.1',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.1',ERRTXT)
       LER=.TRUE.
      ELSEIF(.NOT.(0 .LE. A .AND. A .LT. 1)) THEN
       WRITE(ERRTXT,102) A
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.2',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.2',ERRTXT)
       LER=.TRUE.
      ELSEIF(ABS(NMAX) .GT. 100) THEN
       WRITE(ERRTXT,103) NMAX
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.3',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.3',ERRTXT)
       LER=.TRUE.
      END IF
      IF(LER) RETURN
      EPS=HF*Z10**(-ND)
      NMX=ABS(NMAX)
      IF(NMAX .LE. 0) NMX=1
      DO 5 N = 0,NMX
      RR(N)=0
    5 BA(N)=0
      D=TC(8)*ND+TC(9)
      SX=X
      Q=0
      IF(NMX .GT. 0) THEN
       V=0.5*D/NMX
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 6 I = 2,6
    6   T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       ENDIF
       Q=NMX*T
      ENDIF
#if defined(CERNLIB_DOUBLE)
      F=(HF*X)**A/QGAMMA(1+A)
#endif
#if !defined(CERNLIB_DOUBLE)
      F=(HF*X)**A/DGAMMA(1+A)
#endif
      T=1
      V=TC(10)*D/SX
      IF(LIA) THEN
       F=EXP(X)*F
       V=V-TC(10)
      ENDIF
      IF(LJA .OR. LIA .AND. X .LT. D) THEN
       IF(V .LE. 10) THEN
        T=TC(1)
        DO 7 I = 2,6
    7   T=V*T+TC(I)
       ELSE
        U=LOG(V)-TC(7)
        T=V/(U*(1+(TC(7)-LOG(U))/(1+U)))
       ENDIF
      ENDIF
      NU=1+MAX(Q,TC(11)*SX*T)

      MU=-1
    2 MU=MU+1
      AL=1
      IF(LJA) THEN
       DO 3 N = 1,NU/2
       XN=N
    3  AL=AL*(XN+A)/(XN+1)
       R=0
       S=0
       LEV=.TRUE.
       DO 4 N = 2*(NU/2),1,-1
       XN=N
       XA=XN+A
       R=1/(2*XA/X-R)
       IF(N .LE. NMX) RR(N-1)=R
       IF(LEV) THEN
        AL=AL*(XN+2)/(XA+A)
        S=R*(AL*XA+S)
       ELSE
        S=R*S
       ENDIF
       LEV=.NOT.LEV
    4  CONTINUE
      ELSE
       DO 23 N = 1,NU
       XN=N
   23  AL=AL*(XN+2*A)/(XN+1)
       R=0
       S=0
       DO 24 N = NU,1,-1
       XN=N
       XA=XN+A
       XA2=XA+XA
       R=1/(XA2/X+R)
       IF(N .LE. NMX) RR(N-1)=R
       AL=AL*(XN+1)/(XA+A)
       S=R*(XA2*AL+S)
   24  CONTINUE
      ENDIF
      B(0)=F/(1+S)
      DO 10 N = 0,NMX-1
   10 B(N+1)=RR(N)*B(N)
      DO 11 N = 0,NMX
      IF(ABS(B(N)-BA(N)) .GT. EPS*ABS(B(N))) THEN
       DO 12 M = 0,NMX
   12  BA(M)=B(M)
       NU=NU+5
       IF(MU .LE. 50) GO TO 2
       WRITE(ERRTXT,104) X,A
       IF(LJA) CALL MTLPRT(NAMEJ,'C343.4',ERRTXT)
       IF(LIA) CALL MTLPRT(NAMEI,'C343.4',ERRTXT)
       RETURN
      ENDIF
   11 CONTINUE
      IF(NMAX .LT. 0) THEN
       AL=2/X
       B(1)=AL*A*B(0)+SGN*B(1)
       DO 13 N = 1,-NMAX-1
   13  B(N+1)=AL*(A-N)*B(N)+SGN*B(N-1)
      ENDIF
      RETURN
  101 FORMAT('ILLEGAL ARGUMENT X = ',1P,D15.8)
  102 FORMAT('ILLEGAL ORDER A = ',1P,D15.8)
  103 FORMAT('ILLEGAL NMAX = ',I5)
  104 FORMAT('NO CONVERGENCE FOR X = ',1P,D15.8,' A = ',D15.8,
     1       ' TRY SMALLER ND')
      END
#endif

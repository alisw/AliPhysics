*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:04  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION WWHITM(Z,KA,MU)
#endif
#if !defined(CERNLIB_DOUBLE)
      FUNCTION CWHITM(Z,KA,MU)
#endif
C     Computes the Whittaker function
C      M ka,mu (z) = exp(z/2) * z**(mu+1/2) * 1F1(mu-ka+1/2;1+2*mu;z)
C     using
C      1F1(a;b;z) = ro**(-la-1) * exp(-i*ro) * F la (et,ro) / C la (et)
C     where  F la (et,ro)  is the regular Coulomb Wave Function,
C            a = 1+la+i*et, b = 2*la+2, z = -2*i*ro,
C     and where
C      C la (et) = 2**la * exp(-pi*et/2 + (ln gamma(1+la+i*et) +
C                  ln gamma(1+la-i*et))/2 - ln gamma(2*la+2))
C     is the Gamow factor (cf. Coulomb Wave Functions)
#include "gen/impc64.inc"
#include "gen/defc64.inc"
     +   KA,LA,MU
C     IMPLICIT COMPLEX*16 (A-H,I,K,L,M,O-Z)
#include "gen/def64.inc"
     +   AI,AR,BI,BR,CI,CR,ZI,ZR, R1,DELTA,DLN2,HF,PI,PIH
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
#if defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CWHITM/WWHITM')
#endif
#if !defined(CERNLIB_DOUBLE)
      PARAMETER (NAME = 'CWHITM')
#endif
      DIMENSION F(0:1),G(0:1),FP(0:1),GP(0:1),SG(0:1)

      PARAMETER (DELTA = 1D-13)
      PARAMETER (R1 = 1, HF = R1/2)
      PARAMETER (PI   = 3.14159 26535 89793 238D0, PIH = PI/2)
      PARAMETER (DLN2 = 0.69314 71805 59945 309D0)
      PARAMETER (I = (0,1), IPI = I*PI, IHF= I/2)

#if defined(CERNLIB_DOUBLE)
      GLOGAM(Q)=WLOGAM(Q)
#endif
#if !defined(CERNLIB_DOUBLE)
      GLOGAM(Q)=CLOGAM(Q)
#endif
      A=MU-KA+HF
      B=1+2*MU
      C=MU+KA+HF
      AR=A
      BR=B
      CR=C
      ZR=Z
      AI=-I*A
      BI=-I*B
      CI=-I*C
      ZI=-I*Z
      IF(ZR .LE. 0 .AND. ZI .EQ. 0) THEN
       H=0
       WRITE(ERRTXT,101) ZR
       CALL MTLPRT(NAME,'C328.1',ERRTXT)
      ELSE IF(ABS(BI) .LE. DELTA .AND.
     1        ABS(INT(BR+SIGN(DELTA,BR))+ABS(BR)) .LE .DELTA) THEN
       H=0
       WRITE(ERRTXT,102) (BR-1)/2
       CALL MTLPRT(NAME,'C328.2',ERRTXT)
      ELSE IF(ABS(AI) .LE. DELTA .AND.
     1        ABS(INT(AR+SIGN(DELTA,AR))+ABS(AR)) .LE .DELTA) THEN
       S=1
       Q=1
       N1=AR+SIGN(DELTA,AR)
       DO 1 N = 1,-N1
       Q=Q*((N1+N-1)/(B+(N-1)))*Z/N
    1  S=S+Q
       H=EXP(-HF*Z)*Z**(HF+MU)*S
      ELSE IF(ABS(CI) .LE. DELTA .AND.
     1        ABS(INT(CR+SIGN(DELTA,CR))+ABS(CR)) .LE .DELTA) THEN
       S=1
       Q=1
       N1=CR+SIGN(DELTA,CR)
       DO 2 N = 1,-N1
       Q=-Q*((N1+N-1)/(B+(N-1)))*Z/N
    2  S=S+Q
       H=EXP(HF*Z)*Z**(HF+MU)*S
      ELSE
       ET=I*KA
       RO=IHF*Z
       LA=MU-HF
       Q=GLOGAM(B)-HF*(GLOGAM(A)+GLOGAM(C))+DLN2
       IF(ZR .LE. 0 .AND. ZI .GT. 0) THEN
        ET=-ET
        RO=-RO
        Q=Q+IPI*(MU+HF)
       END IF
#if defined(CERNLIB_DOUBLE)
       CALL WCLBES(RO,ET,LA,0,F,G,FP,GP,SG,-1,4,JFAIL,0)
#endif
#if !defined(CERNLIB_DOUBLE)
       CALL CCLBES(RO,ET,LA,0,F,G,FP,GP,SG,-1,4,JFAIL,0)
#endif
       H=EXP(Q+PIH*(ET-(MU+HF)*I))*F(0)
       IF(JFAIL .NE. 0) THEN
        WRITE(ERRTXT,103) JFAIL
        CALL MTLPRT(NAME,'C328.3',ERRTXT)
       ENDIF
      ENDIF
#if defined(CERNLIB_DOUBLE)
      WWHITM=H
#endif
#if !defined(CERNLIB_DOUBLE)
      CWHITM=H
#endif
      RETURN
  101 FORMAT('ARGUMENT Z =',1P,D15.8,' ON NEGATIVE REAL AXIS')
  102 FORMAT('M KA,MU (Z) UNDEFINED OR INFINITE FOR M = ',F8.1,
     1       ' NEGATIVE HALF-INTEGER')
  103 FORMAT('SUBROUTINE C309 CCLBES/WCLBES RETURNS IFAIL = ',I4)
      END

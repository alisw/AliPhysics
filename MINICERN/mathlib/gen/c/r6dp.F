*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:57  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION C309R6(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
C
C     evaluate the ASYMPTOTIC EXPANSION for the
C            LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION
C
C ***        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)
C
C      that is valid for DREAL(RHO)>0, and best for RHO >> ETA**2, XL,
C      and is derived from the 2F0 expansions for H+ and H-
C      e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)
C      Some lines of this subprogram are for convenience copied from
C           Takemasa, Tamura & Wolter CPC 17 (1979) 351.
C
C     Evaluate to accuracy EPS with at most NMAX terms.
C
C     If the terms start diverging,
C     the corresponding continued fraction is found by RCF
C     & evaluated progressively by Steed's method to obtain convergence.
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION XX(2,NMAX),G(NMAX),C(NMAX)
      DOUBLE PRECISION RE,EPS,T1,T2,T3,AT,ATL,ABSC,FPMAX
      DOUBLE PRECISION HPI

      PARAMETER(HPI  = 1.57079 63267 94896 619D0)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif
      ABSC(W)=ABS(DREAL(W))+ABS(DIMAG(W))
C
      T1=SIN(DREAL(PSI))
      T2=COS(DREAL(PSI))
      ATL=TANH(DIMAG(PSI))

C         GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGN

      COSL=DCMPLX(T2,-T1*ATL)
      TANL=DCMPLX(T1,T2*ATL)/COSL
      RE=0
      XLL1=XL*XL+XL
      ETASQ=ETA*ETA
      SL1=1
      SL=SL1
      SC1=0
      SC=SC1
      TL1=SC
      TL=TL1
      TC1=1-ETA/RHO
      TC=TC1
      FCL=TL+SL*TANL
      G(1)=(TC+SC*TANL)/FCL
      GLAST=G(1)
      ATL=ABSC(GLAST)
      F=GLAST
      D=1
      DF=GLAST
      J=0
      DO 10 N = 2,NMAX
      T1=N-1
      T2=2*T1-1
      T3=T1*T1-T1
      DENOM=2*RHO*T1
      C1=(ETA*T2)/DENOM
      C2=(ETASQ+XLL1-T3)/DENOM
      SL2=C1*SL1-C2*TL1
      TL2=C1*TL1+C2*SL1
      SC2=C1*SC1-C2*TC1-SL2/RHO
      TC2=C1*TC1+C2*SC1-TL2/RHO
      SL=SL+SL2
      TL=TL+TL2
      SC=SC+SC2
      TC=TC+TC2
      SL1=SL2
      TL1=TL2
      SC1=SC2
      TC1=TC2
      FCL=TL+SL*TANL
      IF(ABSC(FCL) .GT. FPMAX .OR. ABSC(FCL) .LT. 1./FPMAX) THEN
       C309R6=G(1)
       FCL=1
       RE=1
       NUSED=0
       RETURN
      END IF
      GSUM=(TC+SC*TANL)/FCL
      G(N)=GSUM-GLAST
      GLAST=GSUM
      AT=ABSC(G(N))
      IF(AT .LT. ABSC(GSUM)*EPS) THEN
       FCL=FCL*COSL
       C309R6=GSUM
       RE=AT/ABSC(GSUM)
       NUSED=N
       RETURN
      END IF
      IF(J .GT. 0 .OR. AT .GT. ATL .OR. N .GE. NMAX-2) J=J+1
      IF(J .EQ. 0) GO TO 10
      CALL C309R7(G,C,J,N,XX,EPS)
      IF(N .LT. 0) THEN
       C309R6=G(1)
       FCL=1
       RE=1
       NUSED=0
       RETURN
      END IF
      DO 60 K = MAX(J,2),N
      D=1/(D*C(K)+1)
      DF=DF*D-DF
      F=F+DF
      IF(ABSC(DF) .LT. ABSC(F)*EPS .OR.
     1   DF .EQ. 0 .AND. F .EQ. 0 .AND. N .GE. 4) THEN
       C309R6=F
       FCL=FCL*COSL
       RE=ABSC(DF)/ABSC(F)
       NUSED=K
       RETURN
      END IF
   60 CONTINUE
      J=N
   10 ATL=AT
      C309R6=F
      FCL=FCL*COSL
      RE=ABSC(DF)/ABSC(F)
      NUSED=-NMAX
      RETURN
      END
#endif

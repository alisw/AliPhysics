#include "isajet/pilot.h"
      SUBROUTINE SIGINT(F,Z,A1S,B1S,A2S,B2S)
C
C          F(N+1) = INT(-Z,Z)(DX X**N/((A1+B1*X)*A2+B2*X)))
C          F(8) = F(9) = 0 (DUMMY VALUES)
C
      DIMENSION F(9)
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION A1,B1,A2,B2,A,B,C,Z,F,A1S,B1S,A2S,B2S
#endif
C
      A1=A1S
      B1=B1S
      A2=A2S
      B2=B2S
      F(8)=0.
      F(9)=0.
C
C          SPECIAL CASE: X**N/(A1*A2)
      IF(B1.EQ.0..AND.B2.EQ.0.) THEN
        F(1)=2.*Z/(A1*A2)
        F(2)=0.
        F(3)=2.*Z**3/(3.*A1*A2)
        F(4)=0.
        F(5)=2.*Z**5/(5.*A1*A2)
        F(6)=0.
        F(7)=2.*Z**7/(7.*A1*A2)
        RETURN
      ENDIF
C
C          SPECIAL CASE: X**N/(A+BX)
      IF(B1.EQ.0..OR.B2.EQ.0.) THEN
        IF(B1.EQ.0.) THEN
          A=A2/B2
          C=1./(A1*B2)
        ELSE
          A=A1/B1
          C=1./(A2*B1)
        ENDIF
        F(1)=LOG((A+Z)/(A-Z))
        F(1)=F(1)*C
        F(2)=-A*LOG((A+Z)/(A-Z))+2.*Z
        F(2)=F(2)*C
        F(3)=A**2*LOG((A+Z)/(A-Z))-2.*A*Z
        F(3)=F(3)*C
        F(4)=-A**3*LOG((A+Z)/(A-Z))+2.*A**2*Z+2.*Z**3/3.
        F(4)=F(4)*C
        F(5)=A**4*LOG((A+Z)/(A-Z))-2.*A**3*Z-2.*A*Z**3/3.
        F(5)=F(5)*C
        F(6)=-A**5*LOG((A+Z)/(A-Z))+2.*A**4*Z+2.*A**2*Z**3/3.+2.*Z**5/5.
        F(6)=F(6)*C
        F(7)=A**6*LOG((A+Z)/(A-Z))-2.*A**5*Z-2.*A**3*Z**3/3.
     $  -2.*A*Z**5/5.
        F(7)=F(7)*C
        RETURN
      ENDIF
C
C          B1 AND B2 NONZERO
      A1=A1/B1
      A2=A2/B2
      C=1./(B1*B2)
C
C          SPECIAL CASE: X**N/(A+B*X)**2
      IF(A1.EQ.A2) THEN
        A=A1
        F(1)=2.*Z/(A**2-Z**2)
        F(1)=F(1)*C
        F(2)=-2.*A*Z/(A**2-Z**2)+LOG((A+Z)/(A-Z))
        F(2)=F(2)*C
        F(3)=(4.*A**2*Z-2.*Z**3)/(A**2-Z**2)-2.*A*LOG((A+Z)/(A-Z))
        F(3)=F(3)*C
        F(4)=(4.*A*Z**3-6.*A**3*Z)/(A**2-Z**2)+3.*A**2*LOG((A+Z)/(A-Z))
        F(4)=F(4)*C
        F(5)=(-16.*A**2*Z**3/3.+8.*A**4*Z-2.*Z**5/3.)/(A**2-Z**2)
     $  -4.*A**3*LOG((A+Z)/(A-Z))
        F(5)=F(5)*C
        F(6)=(4.*A*Z**5/3.+20.*A**3*Z**3/3.-10.*A**5*Z)/(A**2-Z**2)
     $  +5*A**4*LOG((A+Z)/(A-Z))
        F(6)=F(6)*C
        F(7)=(-8.*A**2*Z**5/5.-8.*A**4*Z**3+12.*A**6*Z-2.*Z**7/5.)
     $  /(A**2-Z**2)-6.*A**5*LOG((A+Z)/(A-Z))
        F(7)=F(7)*C
        RETURN
      ENDIF
C
C          GENERAL CASE
      F(1)=(-LOG((A1+Z)/(A1-Z))+LOG((A2+Z)/(A2-Z)))/(A1-A2)
      F(1)=F(1)*C
      F(2)=(A1*LOG((A1+Z)/(A1-Z))-A2*LOG((A2+Z)/(A2-Z)))/(A1-A2)
      F(2)=F(2)*C
      F(3)=(-A1**2*LOG((A1+Z)/(A1-Z))+A2**2*LOG((A2+Z)/(A2-Z)))/(A1-A2)
     $+2.*Z
      F(3)=F(3)*C
      F(4)=(A1**3*LOG((A1+Z)/(A1-Z))-A2**3*LOG((A2+Z)/(A2-Z)))/(A1-A2)
     $+2.*Z*(-A1-A2)
      F(4)=F(4)*C
      F(5)=(-A1**4*LOG((A1+Z)/(A1-Z))+A2**4*LOG((A2+Z)/(A2-Z)))/(A1-A2)
     $+2.*Z*(A1*A2+A1**2+A2**2)+2.*Z**3/3.
      F(5)=F(5)*C
      F(6)=(A1**5*LOG((A1+Z)/(A1-Z))-A2**5*LOG((A2+Z)/(A2-Z)))/(A1-A2)
     $+2.*Z*(-A1*A2**2-A1**2*A2-A1**3-A2**3)+2.*Z**3/3.*(-A1-A2)
      F(6)=F(6)*C
      F(7)=(-A1**6*LOG((A1+Z)/(A1-Z))+A2**6*LOG((A2+Z)/(A2-Z)))/(A1-A2)
     $+2.*Z*(A1*A2**3+A1**2*A2**2+A1**3*A2+A1**4+A2**4)
     $+2.*Z**3/3.*(A1*A2+A1**2+A2**2)+2.*Z**5/5.
      F(7)=F(7)*C
      RETURN
      END

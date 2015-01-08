C--------------------------------------------------------------------------
C
C Environment:
C      This software is part of the EvtGen package developed jointly
C      for the BaBar and CLEO collaborations.  If you use all or part
C      of it, please give an appropriate acknowledgement.
C
C Copyright Information: See EvtGen/COPYRIGHT
C      Copyright (C) 1998      Caltech, UCSB
C
C Module: EvtLi2Spence.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
C -----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION LI2SPENCE(X)
C *************************************************************
C REAL SPENCE FUNCTION OBTAINED FROM J. VAN DER BIJ
C *************************************************************
       IMPLICIT REAL*8 (A-Z)
      COMPLEX*16 Y,Z,Z2,SP
      REAL*8 XR,XI,R2
      REAL*8 A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,ZETA2
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,ZETA2/
     1 -0.250000000000000D0,
     2 -0.111111111111111D0,
     3 -0.010000000000000D0,
     4 -0.017006802721088D0,
     5 -0.019444444444444D0,
     6 -0.020661157024793D0,
     7 -0.021417300648069D0,
     8 -0.021948866377231D0,
     9 -0.022349233811171D0,
     1 -0.022663689135191D0,
     2  1.644934066848226D0/

      LI2SPENCE=0
      XR=X
      XI=0.0D0
      R2=XR*XR+XI*XI
      IF(R2.GT.1.D0.AND.(XR/R2).GT.0.5D0)THEN
        Y=(X-1.D0)/X
        Z=-LOG(1.D0-Y)
        Z2=Z*Z
        SP=Z*(1.D0+A1*Z*(1.D0+A2*Z*(1.D0+A3*Z2*(1.D0+A4*Z2*
     1 (1.D0+A5*Z2*(1.D0+A6*Z2*(1.D0+A7*Z2*(1.D0+A8*Z2*(1.D0+A9*Z2*
     2 (1.D0+A10*Z2))))))))))
     3 +ZETA2-LOG(X)*LOG(1.D0-X)+0.5D0*LOG(X)**2
        LI=DREAL(SP)
        RETURN
            ELSEIF(R2.GT.1.D0.AND.(XR/R2).LE.0.5D0)THEN
        Y=1.D0/X
        Z=-LOG(1.D0-Y)
        Z2=Z*Z
        SP=-Z*(1.D0+A1*Z*(1.D0+A2*Z*(1.D0+A3*Z2*(1.D0+A4*Z2*
     1 (1.D0+A5*Z2*(1.D0+A6*Z2*(1.D0+A7*Z2*(1.D0+A8*Z2*(1.D0+A9*Z2*
     2 (1.D0+A10*Z2))))))))))
     3 -ZETA2-0.5D0*LOG(-X)**2
        LI=DREAL(SP)
        RETURN
      ELSEIF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        SP=ZETA2
        LI=DREAL(SP)
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        Z=-LOG(1.D0-Y)
        Z2=Z*Z
        SP=-Z*(1.D0+A1*Z*(1.D0+A2*Z*(1.D0+A3*Z2*(1.D0+A4*Z2*
     1 (1.D0+A5*Z2*(1.D0+A6*Z2*(1.D0+A7*Z2*(1.D0+A8*Z2*(1.D0+A9*Z2*
     2 (1.D0+A10*Z2))))))))))
     3 +ZETA2-LOG(X)*LOG(1.D0-X)
       LI=DREAL(SP)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        Z=-LOG(1.D0-Y)
        Z2=Z*Z
        SP=Z*(1.D0+A1*Z*(1.D0+A2*Z*(1.D0+A3*Z2*(1.D0+A4*Z2*
     1 (1.D0+A5*Z2*(1.D0+A6*Z2*(1.D0+A7*Z2*(1.D0+A8*Z2*(1.D0+A9*Z2*
     2 (1.D0+A10*Z2))))))))))
        LI2SPENCE=DREAL(SP)
        RETURN
      ENDIF
      END
C ----------------------------------------------------------------------
 

*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      FUNCTION SHFPI(X)
c	=================

c	Pion Pt parametrization: mt at low pt, CDF at high pt
c	[ANL-HEP-PR 88-32]

      P0=1.3
      XN=8.28
      XLIM=0.5
      T=.160
      AM=ULMASS(111)
      AM2=AM**2
      IF (X.LT.XLIM) THEN
        SHFPI=X*EXP(-SQRT(X**2+AM2)/T)
      ELSE
        SHFPI=X*EXP(-SQRT(XLIM**2+AM2)/T)*((P0+XLIM)**XN)/((P0+X)**XN)
      ENDIF

      RETURN
      END

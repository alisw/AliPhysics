*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      FUNCTION SHMTSC(AM,X)
c	=====================

c	Mt scaling for a particle of mass AM
c	[E.g.: V.Hedberg, LUNFD6/(NFFL-7073)/1987;
c       M.Bourquin and J.M.Gaillard, Nucl. Phys. B114 (1976),334]

      B = 2.	
      XM = 12.3

      AMPI0 = ULMASS(111)
      TMPI = SQRT(AMPI0**2+X**2)
      TMPART = SQRT(AM**2+X**2)
      SHMTSC = ((TMPI+B)/(TMPART+B))**XM * SHFPI(X)

      RETURN
      END

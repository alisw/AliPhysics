*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      FUNCTION SHFPRO(X)
c	==================

c	Proton pt distributions from Tevatron
c	[T.Alexopoulos et al.: Phys. Rev. Lett. 64 (1990), 991]

      A = 2.8

      SHFPRO = X*EXP(-A*X)

      RETURN
      END

*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      FUNCTION SHFKAO(X)
c	==================

c	Kaon pt distributions from Tevatron
c	[T.Alexopoulos et al.: Phys. Rev. Lett. 64 (1990), 991]

      A = 3.69

      SHFKAO = X*EXP(-A*X)

      RETURN
      END

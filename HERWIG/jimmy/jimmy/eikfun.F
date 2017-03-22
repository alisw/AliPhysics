*CMZ :          25/04/95  15.01.42  by  Unknown
*-- Author : Jeff Forshaw
      DOUBLE PRECISION FUNCTION EIKFUN(NDIM,V)
C ----------------------------------------------------------------------------
C  Purpose: Used when eikonalising the total gamma-p xsec.
C ----------------------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"


      INTEGER NDIM

      DOUBLE PRECISION XSECN,V(NDIM),B2MAX,B2,FAC1
      DOUBLE PRECISION JMAREA, AREA
      COMMON/BLOCK1/XSECN

      B2MAX=5.D2
      FAC1=B2MAX
      B2=FAC1*V(1)
      AREA = JMAREA(B2)
      EIKFUN=PIFAC/PHAD*(1.D0-DEXP(-AREA*XSECN*PHAD))*FAC1

      RETURN
      END

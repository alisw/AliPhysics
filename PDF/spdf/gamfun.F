*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:27  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
C     Gamma function : See Abramowitz, page 257, form. 6.4.40
C
#if defined(CERNLIB_SINGLE)
      FUNCTION GAMFUN(Y)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION FUNCTION GAMFUN(Y)
#endif
#include "pdf/impdp.inc"
#include "pdf/expdp.inc"
     +     Y,R,S,T,AFSPL,X,
     +     COEF(10),PI,ZEROD,HALFD,ONED,TWOD,TEND
C
      DATA COEF/8.3333333333333334D-02,-2.7777777777777778D-03,
     .          7.9365079365079365D-04,-5.9523809523809524D-04,
     .          8.4175084175084175D-04,-1.9175269175269175D-03,
     .          6.4102564102564103D-03,-2.9550653594771242D-02,
     .          0.1796443723688306D0  ,-0.6962161084529506D0  /
      DATA PI/  3.141592653589793D0/
      DATA ZEROD/0.D0/,HALFD/0.5D0/,ONED/1.D0/,TWOD/2.D0/,TEND/10.D0/
C
      X=Y
      AFSPL=ONED
      N=INT(TEND-Y)
      DO 10 I=0,N
        AFSPL=AFSPL*X
        X=X+ONED
10    CONTINUE
      R=(X-HALFD)* LOG(X)-X+HALFD* LOG(TWOD*PI)
      S=X
      T=ZEROD
      DO 20 I=1,10
        T=T+COEF(I)/S
        S=S*X**2
20    CONTINUE
      GAMFUN = EXP(R+T)/AFSPL
      END

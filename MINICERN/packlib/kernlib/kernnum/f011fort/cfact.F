*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          CFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
          INTEGER             IR(*),    IPAIRF
          COMPLEX             A(IDIM,*),DET,      ZERO,     ONE,X,Y,TF
          REAL                G1,       G2
          REAL                PIVOTF,   P,        Q,        SIZEF,  T
          CHARACTER*6         HNAME
#if !defined(CERNLIB_NUMCDBLE)
          COMPLEX             S11, S12, DOTF
          DOTF(X,Y,S11)  =  X*Y + S11
#endif
#if defined(CERNLIB_NUMCDBLE)
          COMPLEX*16          S11, S12, DC, DOTF
          DC(X)        =  DCMPLX(DBLE(REAL(X)),DBLE(AIMAG(X)))
          DOTF(X,Y,S11)  =  DC(X)*DC(Y) + S11
#endif
          IPAIRF(J,K)  =  J*2**12 + K
          PIVOTF(X)    =  MAX(ABS(REAL(X)),ABS(AIMAG(X)))
          SIZEF(X)     =  MAX(ABS(REAL(X)),ABS(AIMAG(X)))
#if defined(CERNLIB_NUME2465)
          DATA      G1, G2              /  1.E-1232, 1.E1232  /
#endif
#if defined(CERNLIB_NUME293)
          DATA      G1, G2              /  1.E-146, 1.E146  /
#endif
#if defined(CERNLIB_NUME75)
          DATA      G1, G2              /  1.E-37,  1.E37  /
#endif
#if defined(CERNLIB_NUME38)
          DATA      G1, G2              /  1.E-19,  1.E19  /
#endif
#if defined(CERNLIB_NUME999)
          DATA      ?????  G1, G2 NOT DEFINED  ?????
#endif
          DATA      HNAME               /  ' CFACT'  /
          DATA      ZERO, ONE           /  (0., 0.), (1., 0.)  /
          DATA      NORMAL, IMPOSS      /  0, -1  /
          DATA      JRANGE, JOVER, JUNDER  /  0, +1, -1  /
#include "fact.inc"
          RETURN
          END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          DFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
          INTEGER             IR(*),    IPAIRF
          DOUBLE PRECISION    A(IDIM,*),DET,      ZERO,     ONE,X,Y,TF
          REAL                G1,       G2
          REAL                PIVOTF,   P,        Q,        SIZEF,  T
          DOUBLE PRECISION    S11, S12, DOTF
          CHARACTER*6         HNAME
          IPAIRF(J,K)  =  J*2**12 + K
          PIVOTF(X)    =  ABS(SNGL(X))
          SIZEF(X)     =  ABS(SNGL(X))
          DOTF(X,Y,S11)  =  X * Y + S11
#if defined(CERNLIB_NUME2465)
          DATA      G1, G2              /  1.E-1232, 1.E1232  /
#endif
#if defined(CERNLIB_NUME279)
          DATA      G1, G2              /  1.E-139, 1.E139  /
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
          DATA      HNAME               /  ' DFACT'  /
          DATA      ZERO, ONE           /  0.D0, 1.D0  /
          DATA      NORMAL, IMPOSS      /  0, -1  /
          DATA      JRANGE, JOVER, JUNDER  /  0, +1, -1  /
#include "fact.inc"
          RETURN
          END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          CFINV(N,A,IDIM,IR)
          INTEGER             IR(*)
          COMPLEX             A(IDIM,*),     X, Y, TI
          CHARACTER*6         HNAME
#if !defined(CERNLIB_NUMCDBLE)
          COMPLEX             ZERO, S31, S32, S33, S34, DOTF
          DOTF(X,Y,S31)  =  X*Y + S31
          DATA      ZERO      /  (0.,0.)  /
#endif
#if defined(CERNLIB_NUMCDBLE)
          COMPLEX*16          ZERO, S31, S32, S33, S34, DC, DOTF
          DC(X)          =  DCMPLX(DBLE(REAL(X)),DBLE(AIMAG(X)))
          DOTF(X,Y,S31)  =  DC(X)*DC(Y) + S31
          DATA      ZERO      /  (0.D0,0.D0)  /
#endif
          DATA      HNAME               /  ' CFINV'  /
#include "finv.inc"
          RETURN
          END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:48  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION G116F1(X)

      COMMON /G116C1/ H(7),T0,T1,T,OMEGA,A(155),B(155),X0

      G116F1=H(1)+H(2)*LOG(H(3)*X)-H(4)*X
      RETURN

      ENTRY G116F2(X)

      G116F2=H(5)-X+H(6)*(LOG(ABS(X))+REXPIN(X))-H(7)*EXP(-X)
      RETURN
      END

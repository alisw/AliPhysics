*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UCOCOP(A,B,IDO,IW,NA,NB)
C
C CERN PROGLIB# V302    UCOCOP          .VERSION KERNFOR  4.40  940929
C ORIG. 01/01/70 F.BRUYANT
C
#if defined(CERNLIB_QINTCOPY)
      INTEGER   A(*),B(*)
#endif
#if !defined(CERNLIB_QINTCOPY)
      DIMENSION A(*),B(*)
#endif
C
C
      IF (IDO.EQ.0) RETURN
      IF (IW.EQ.0) RETURN
      K=0
      L=0
      KK=NA-IW
      LL=NB-IW
      DO 20 I=1,IDO
      DO 10 J=1,IW
      K=K+1
      L=L+1
   10 B(L)=A(K)
      K=K+KK
   20 L=L+LL
      RETURN
      END

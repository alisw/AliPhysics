*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UDICOP(A,B,IDO,IW,NA,NB)
C
C CERN PROGLIB# V302    UDICOP          .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71
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
      K=IW+(IDO-1)*NA
      L=IW+(IDO-1)*NB
      KK=NA-IW
      LL=NB-IW
         DO 40 I=1,IDO
         DO 30 J=1,IW
      B(L)=A(K)
      K=K-1
   30 L=L-1
      K=K-KK
   40 L=L-LL
      RETURN
      END

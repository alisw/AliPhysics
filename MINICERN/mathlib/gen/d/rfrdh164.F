*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:22  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DFRDH1(F,G,M,TV,NGTV,WS,IDIM,N)
C
#include "gen/imp64.inc"
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='DFRDH1')
#endif
#if !defined(CERNLIB_DOUBLE)
      SUBROUTINE RFRDH1(F,G,M,TV,NGTV,WS,IDIM,N)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='RFRDH1')
#endif
C
      DIMENSION TV(0:*),NGTV(*),WS(IDIM,*)
      CHARACTER*80 ERRTXT

      N0=0
      DO 1 K = 1,M
    1 N0=N0+NGTV(K)
      N1=1
      DO 2 K = 1,M
#if defined(CERNLIB_DOUBLE)
      CALL DGSET(TV(K-1),TV(K),NGTV(K),WS(N1,N0+3),WS(N1,N0+4))
#endif
#if !defined(CERNLIB_DOUBLE)
      CALL RGSET(TV(K-1),TV(K),NGTV(K),WS(N1,N0+3),WS(N1,N0+4))
#endif
    2 N1=N1+NGTV(K)

      DO 3 I = 1,N0
      DO 4 J = 1,N0
    4 WS(I,J)=-WS(J,N0+4)*G(WS(I,N0+3),WS(J,N0+3))
      WS(I,I)=1+WS(I,I)
    3 WS(I,N0+1)=F(WS(I,N0+3))
#if defined(CERNLIB_DOUBLE)
      CALL DEQN(N0,WS,IDIM,WS(1,N0+2),IFAIL,1,WS(1,N0+1))
#endif
#if !defined(CERNLIB_DOUBLE)
      CALL REQN(N0,WS,IDIM,WS(1,N0+2),IFAIL,1,WS(1,N0+1))
#endif
      IF(IFAIL .EQ. -1) THEN
       CALL MTLPRT(NAME,'D601.1','MATRIX SINGULAR')
      ENDIF
      N=N0
      RETURN
      END

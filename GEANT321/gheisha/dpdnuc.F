*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:06  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      FUNCTION DPDNUC(A,B,C)
C
C CALLED BY : PHPNUC
C ORIGIN    : H.FESEFELDT
C
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
C
C     DPDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
      A2 = A*A
      B2 = B*B
      C2 = C*C
      IF(A2) 21,21,61
   61 CONTINUE
      ARG=A2+(B2-C2)**2/A2-2.0*(B2+C2)
      IF (ARG) 21,21,31
   21 DPDNUC=0.0
      GOTO 41
   31 CONTINUE
      DPDNUC = 0.5*SQRT(A2 + (B2-C2)**2/A2 - 2.0*(B2+C2))
   41 CONTINUE
      RETURN
      END

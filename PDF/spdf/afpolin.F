*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:33  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:00  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE AFPOLIN(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10)
      DIMENSION XA(NMAX),YA(NMAX),C(NMAX),D(NMAX)
      Y=0.
      IF(N.GT.NMAX) RETURN
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
C         IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

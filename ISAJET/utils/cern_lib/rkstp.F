*
* $Id$
*
* $Log$
* Revision 1.9  2000/07/25 14:53:05  mclareni
* Version 7.51 from author
*
*
*#include "sys/CERNLIB_machine.h"
#include "isajet/pilot.h"
#if defined(CERNLIB_NOCERN)
C-----------------------------------------------------------------------
      SUBROUTINE RKSTP(N,H,X,Y,SUB,W)
C-----------------------------------------------------------------------
C
C     From CERN Program Library, routine D209, with error message for
C     N.LT.1 replaced by STOP 99 to eliminate Kernlib error routine.
C
      DIMENSION Y(N),W(N,3)
      LOGICAL MFLAG,RFLAG
      EXTERNAL SUB
C
C     ******************************************************************
C
C     THIS SUBROUTINE REPLACES X BY X+H AND ADVANCES THE SOLUTION OF THE
C     SYSTEM OF DIFFERENTIAL EQUATIONS DY/DX=F(X,Y) FROM Y(X) TO Y(X+H)
C     USING A FIFTH-ORDER RUNGE-KUTTA METHOD.
C
C     SUB IS THE NAME OF A SUBROUTINE SUB(X,Y,F) WHICH SETS THE VECTOR F
C     TO THE DERIVATIVE AT X OF THE VECTOR Y.
C
C     W IS A WORKING-SPACE ARRAY, TREATED AS CONSISTING OF THREE CONSEC-
C     UTIVE WORKING VECTORS OF LENGTH N.
C
C     ******************************************************************
C
C  START.
      IF (N.LT.1) STOP 99
      NLOCAL=N
      HLOCAL=H
      H2=0.5*HLOCAL
      H6=HLOCAL/6.
      XH=X+HLOCAL
      XH2=X+H2
      CALL SUB(X,Y,W(1,1))
      DO 1 J=1,NLOCAL
         W(J,2)=Y(J)+H2*W(J,1)
    1 CONTINUE
      CALL SUB(XH2,W(1,2),W(1,3))
      DO 2 J=1,NLOCAL
         W(J,1)=W(J,1)+2.*W(J,3)
         W(J,2)=Y(J)+H2*W(J,3)
    2 CONTINUE
      CALL SUB(XH2,W(1,2),W(1,3))
      DO 3 J=1,NLOCAL
         W(J,1)=W(J,1)+2.*W(J,3)
         W(J,2)=Y(J)+HLOCAL*W(J,3)
    3 CONTINUE
      CALL SUB(XH,W(1,2),W(1,3))
      DO 4 J=1,NLOCAL
         Y(J)=Y(J)+H6*(W(J,1)+W(J,3))
    4 CONTINUE
      X=XH
      RETURN
      END
#endif

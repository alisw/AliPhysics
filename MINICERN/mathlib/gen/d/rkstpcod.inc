*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:16  mclareni
* Mathlib gen
*
*
*
* rkstpcod.inc
*
      DIMENSION Y(*),W(N,*)

      PARAMETER (Z1 = 1, HF = Z1/2, HS = Z1/6)

      IF(N .LE. 0) RETURN
      H2=HF*H
      H6=HS*H
      XH=X+H
      XH2=X+H2
      CALL SUB(X,Y,W(1,1))
      DO 1 J=1,N
    1 W(J,2)=Y(J)+H2*W(J,1)
      CALL SUB(XH2,W(1,2),W(1,3))
      DO 2 J=1,N
      W(J,1)=W(J,1)+2*W(J,3)
    2 W(J,2)=Y(J)+H2*W(J,3)
      CALL SUB(XH2,W(1,2),W(1,3))
      DO 3 J=1,N
      W(J,1)=W(J,1)+2*W(J,3)
    3 W(J,2)=Y(J)+H*W(J,3)
      CALL SUB(XH,W(1,2),W(1,3))
      DO 4 J=1,N
    4 Y(J)=Y(J)+H6*(W(J,1)+W(J,3))
      X=XH
      RETURN
      END

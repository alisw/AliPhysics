*
* $Id$
*
* $Log$
* Revision 1.2  1997/12/15 16:18:42  mclareni
* Changes for the Portland Group f77 compiler inside cpp define CERNLIB_QFPGF77
*
* Revision 1.1  1996/04/16 15:57:01  mclareni
* The name of cfstft was mistyped and also becomes D706, not 705
*
*
*
      SUBROUTINE CFSTFT(MS,A)
 
      COMPLEX A(0:*),U,W,T
 
      IF(MS .EQ. 0) GO TO 3
      M=ABS(MS)
      N=2**M
      J=0
      DO 7 I = 0,N-2
      IF(I .LT. J) THEN
       T=A(J)
       A(J)=A(I)
       A(I)=T
      ENDIF
      K=N/2
    6 IF(K .LE. J) THEN
       J=J-K
       K=K/2
       GO TO 6
      ENDIF
    7 J=J+K
      DO 8 I = 0,N-1,2
      T=A(I+1)
      A(I+1)=A(I)-T
      A(I)=A(I)+T
    8 CONTINUE
      C=0
      S=SIGN(1,MS)
      LE=2
      DO 2 L = 2,M
      W=CMPLX(C,S)
      U=W
      C=SQRT(0.5*C+0.5)
#if defined(CERNLIB_QFPGF77)
      S=AIMAG(W)/(C+C)
#else
      S=IMAG(W)/(C+C)
#endif
      LE1=LE
      LE=LE1+LE1
      DO 9 I = 0,N-1,LE
      T=A(I+LE1)
      A(I+LE1)=A(I)-T
      A(I)=A(I)+T
    9 CONTINUE
      DO 2 J = 2,LE1
      DO 1 I = J-1,N-1,LE
      T=A(I+LE1)*U
      A(I+LE1)=A(I)-T
      A(I)=A(I)+T
    1 CONTINUE
      U=U*W
    2 CONTINUE
    3 RETURN
      END

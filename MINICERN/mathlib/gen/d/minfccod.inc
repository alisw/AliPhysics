*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:20  mclareni
* Mathlib gen
*
*
*
* minfccod.inc
*
      PARAMETER (W5 = 2.23606 79774 99790D0)
      PARAMETER (HV = (3-W5)/2, HW = (W5-1)/2, R1 = 1, HF = R1/2)

      N=-1
      IF(A .NE. B) N=NINT(2.08*LOG(ABS((A-B)/EPS)))
      C=A
      D=B
      IF(A .GT. B) THEN
       C=B
       D=A
      ENDIF
      LLT=.TRUE.
      LGE=.TRUE.

    1 H=D-C
      IF(N .LT. 0) THEN
       X=HF*(C+D)
       Y=F(X)
       LLM=ABS(X-A) .GT. DELTA .AND. ABS(X-B) .GT. DELTA
       RETURN
      ENDIF
      IF(LLT) THEN
       V=C+HV*H
       FV=F(V)
      ENDIF
      IF(LGE) THEN
       W=C+HW*H
       FW=F(W)
      ENDIF
      IF(FV .LT. FW) THEN
       LLT=.TRUE.
       LGE=.FALSE.
       D=W
       W=V
       FW=FV
      ELSE
       LLT=.FALSE.
       LGE=.TRUE.
       C=V
       V=W
       FV=FW
      ENDIF
      N=N-1
      GO TO 1
      END

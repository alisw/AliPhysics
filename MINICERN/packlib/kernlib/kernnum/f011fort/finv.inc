*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
*
* finv.inc
*
          IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 310
             CALL TMPRNT(HNAME,N,IDIM,0)
             RETURN
 310      IF(N .EQ. 1)  RETURN
          A(2,1)  =  -A(2,2) * DOTF(A(1,1),A(2,1),ZERO)
          A(1,2)  =  -A(1,2)
          IF(N .EQ. 2)  GOTO 330
          DO 314    I  =  3, N
             IM2  =  I-2
             DO 312 J  =  1, IM2
                S31  =  ZERO
                S32  =  A(J,I)
                DO 311  K  =  J, IM2
                   S31  =  DOTF(A(K,J),A(I,K),S31)
                   S32  =  DOTF(A(J,K+1),A(K+1,I),S32)
 311               CONTINUE
                A(I,J)  =  -A(I,I) * DOTF(A(I-1,J),A(I,I-1),S31)
                A(J,I)  =  -S32
 312            CONTINUE
             A(I,I-1)  =  -A(I,I) * DOTF(A(I-1,I-1),A(I,I-1),ZERO)
             A(I-1,I)  =  -A(I-1,I)
 314         CONTINUE
 330      NM1  =  N-1
          DO 335   I  =  1, NM1
             NMI  =  N-I
             DO 332   J  =  1, I
                S33  =  A(I,J)
                DO 331   K  =  1, NMI
                   S33  =  DOTF(A(I+K,J),A(I,I+K),S33)
 331               CONTINUE
                A(I,J)  =  S33
 332            CONTINUE
             DO 334   J  =  1, NMI
                S34  =  ZERO
                DO 333   K  =  J, NMI
                   S34  =  DOTF(A(I+K,I+J),A(I,I+K),S34)
 333               CONTINUE
                A(I,I+J)  =  S34
 334            CONTINUE
 335         CONTINUE
          NXCH  =  IR(N)
          IF(NXCH .EQ. 0)  RETURN
            DO 342 M  =  1, NXCH
             K   =  NXCH - M+1
             IJ  =  IR(K)
             I   =  IJ / 4096
             J   =  MOD(IJ,4096)
             DO 341  K  =  1, N
                TI      =  A(K,I)
                A(K,I)  =  A(K,J)
                A(K,J)  =  TI
 341            CONTINUE
 342         CONTINUE

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
*
* fact.inc
*
          IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 110
             CALL TMPRNT(HNAME,N,IDIM,0)
             RETURN
 110      IFAIL  =  NORMAL
          JFAIL  =  JRANGE
          NXCH   =  0
          DET    =  ONE
          DO 144    J  =  1, N
 120         K  =  J
             P  =  PIVOTF(A(J,J))
             IF(J .EQ. N)  GOTO 122
             JP1  =  J+1
             DO 121    I  =  JP1, N
                Q  =  PIVOTF(A(I,J))
                IF(Q .LE. P)  GOTO 121
                   K  =  I
                   P  =  Q
 121            CONTINUE
             IF(K .NE. J)  GOTO 123
 122         IF(P .GT. 0.)  GOTO 130
                DET    =  ZERO
                IFAIL  =  IMPOSS
                JFAIL  =  JRANGE
                RETURN
 123         DO 124    L  =  1, N
                TF      =  A(J,L)
                A(J,L)  =  A(K,L)
                A(K,L)  =  TF
 124            CONTINUE
             NXCH      =  NXCH + 1
             IR(NXCH)  =  IPAIRF(J,K)
 130         DET     =  DET * A(J,J)
             A(J,J)  =  ONE / A(J,J)
             T  =  SIZEF(DET)
             IF(T .LT. G1)  THEN
                DET    =  ZERO
                IF(JFAIL .EQ. JRANGE)  JFAIL  =  JUNDER
             ELSEIF(T .GT. G2)  THEN
                DET    =  ONE
                IF(JFAIL .EQ. JRANGE)  JFAIL  =  JOVER
             ENDIF
             IF(J .EQ. N)  GOTO 144
             JM1  =  J-1
             JP1  =  J+1
             DO 143   K  =  JP1, N
                S11  =  -A(J,K)
                S12  =  -A(K,J+1)
                IF(J .EQ. 1)  GOTO 142
                DO 141  I  =  1, JM1
                   S11  =  DOTF(A(I,K),A(J,I),S11)
                   S12  =  DOTF(A(I,J+1),A(K,I),S12)
 141               CONTINUE
 142            A(J,K)    =  -S11 * A(J,J)
                A(K,J+1)  =  -DOTF(A(J,J+1),A(K,J),S12)
 143            CONTINUE
 144         CONTINUE
 150      IF(MOD(NXCH,2) .NE. 0)  DET  =  -DET
          IF(JFAIL .NE. JRANGE)   DET  =  ZERO
          IR(N)  =  NXCH

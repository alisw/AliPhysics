*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:41  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE BITPOS(I,N,K,M)
C
C CERN PROGLIB# M508    BITPOS          .VERSION KERNFOR  4.16  870601
C ORIG.  OCT 81,  M.METCALF, CERN/DD
C
C     TO INDICATE WHICH BITS IN A SERIES OF WORDS ARE SET.
C     BITS WITHIN A WORD ARE NUMBERED RIGHT-TO-LEFT, STARTING AT 0.
C
C     I=INPUT WORDS
C     N=NO. OF BITS TO BE TESTED
C     K=POSITION ARRAY
C     M=NO. OF SET BITS IN THE FIRST N POSITIONS
C
#include "kerngen/wordsize.inc"
      PARAMETER   (LBIT = NBITPW)
C
      INTEGER I(*),K(*)
      LOGICAL BTEST
C
C     INITIALIZE
      NWORD=(N-1)/LBIT+1
      M=0
      JND=1
      NBIT = LBIT
C
C     UNPACK EACH WORD
      DO 1 MM=1,NWORD
         NSET = 0
         IMM = I(MM)
         IF(MM.EQ.NWORD) NBIT=N-(NWORD-1)*LBIT
         JU = JND-1
C
C     LOCATE SET BITS
         DO 4 MMU = 1,NBIT
            IF (BTEST(IMM,0)) THEN
               JU = JU+1
               K(JU)=MMU-1
               NSET = NSET+1
            ENDIF
         IMM = ISHFT(IMM,-1)
    4    CONTINUE
    9    M = JU
         IF(MM.EQ.1) GO TO 2
C
C     BIAS WORDS OTHER THAN THE FIRST
         KND=JND+NSET-1
         IBIAS=(MM-1)*LBIT
         DO 3 NN=JND,KND
    3    K(NN)=K(NN)+IBIAS
    2    JND=JND+NSET
    1 CONTINUE
      END

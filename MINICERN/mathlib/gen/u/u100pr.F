*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:46  mclareni
* Mathlib gen
*
*
*FCA :          Fri Mar 26 17:27:50 CET 1999 by  Federico Carminati
*               removed vzero
#include "gen/pilot.h"
      SUBROUTINE U100PR (I,K,M)
      DIMENSION IPRIME(40),K(40)
      DATA IPRIME/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
     + 71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
     + 157,163,167,173/
      II=I
      M=1
      IF(I.LT.1)  GO TO 6
      DO KZERO=1,40
         K(KZERO)=0
      ENDDO
*      CALL VZERO (K, 40)
      IF(I.EQ.1)  RETURN
      DO 2 KL=1,40
    1 IF(MOD(I,IPRIME(KL)).NE.0)   GO TO 2
      K(KL)=K(KL)+1
      I=I/IPRIME(KL)
      IF(I.EQ.1)  GO TO 3
      GO TO 1
    2 CONTINUE
      IF(I.NE.1)    WRITE(6,4)    I
    4 FORMAT(' PRIME NO. LARGER THAN THE LIST,THE FACTOR LEFT=',I15)
      M=I
    3 I=II
      RETURN
    6 WRITE(6,7)
    7 FORMAT(' PLEASE OBSERVE THE CONVENTION OF YOUR ASKING NO. BEING',
     F'LARGER THAN +0')
      RETURN
      END

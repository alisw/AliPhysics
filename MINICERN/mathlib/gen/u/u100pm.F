*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:45  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE U100PM(J,IN,ID)
      DIMENSION  J(40),IPRIME(40)
      DATA IPRIME/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
     + 71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
     + 157,163,167,173/
      IN=1
      ID=1
      DO 1 LL=1,40
      IF(J(LL).GT.0)    IN=IN*IPRIME(LL)**J(LL)
    1 IF(J(LL).LT.0)    ID=ID*IPRIME(LL)**(-J(LL))
      RETURN
      END

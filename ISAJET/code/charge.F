#include "isajet/pilot.h"
      FUNCTION CHARGE(ID)
C
C          COMPUTE CHARGE OF PARTICLE WITH IDENT CODE ID
C          ICHRG MUST BE DIMENSIONED NQLEP+13
C
#include "isajet/itapes.inc"
      DIMENSION ICHRG(75),IFL(3)
C          3 * charge
      DATA ICHRG/0
     $,2,-1,-1,2,-1,2,-1,2,0,0, 0,-3,0,-3,0,-3,0,-3,0,0,0
     $,2,-1,-1,2,-1,2,-1,2,0,0, 0,-3,0,-3,0,-3,0,-3,3,0
     $,2,-1,-1,2,-1,2,-1,2,3,0, 0,-3,0,-3,0,-3,0,-3,3,0
     $,3,0,0,0,0,0,3,3,6,6,0,0,0/
C
      IDABS=IABS(ID)
      CALL FLAVOR(ID,IFL(1),IFL(2),IFL(3),JSPIN,INDEX)
      IF(IDABS.LT.100) GO TO 200
C
      ISUM=0
      DO 100 I=1,3
        ISUM=ISUM+ICHRG(IABS(IFL(I))+1)*ISIGN(1,IFL(I))
  100 CONTINUE
      CHARGE=ISUM/3.
      RETURN
C
200   CHARGE=ICHRG(INDEX+1)*ISIGN(1,ID)
      CHARGE=CHARGE/3.
      RETURN
      END

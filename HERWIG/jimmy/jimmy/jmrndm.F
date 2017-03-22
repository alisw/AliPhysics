      DOUBLE PRECISION FUNCTION JMRNDM(WOT,SEED)
C ------------------------------------------------------------------------
C  Random number generator interface.
C  1) SET  2) GENERATE UNIFORM RANDOM NUMBER 3) RETURN CURRENT SEEDS
C ------------------------------------------------------------------------

      IMPLICIT NONE

#include "jimmy.inc"


      INTEGER WOT, SEED(2)
      DOUBLE PRECISION HWRGEN, HWRGET, HWRSET

      IF (WOT.EQ.1) THEN

C       Initialise
        JMRNDM = HWRSET(SEED)

      ELSE IF (WOT.EQ.2) THEN

C       Generate a random number
        JMRNDM = HWRGEN(SEED)

      ELSE IF (WOT.EQ.3) THEN

C       Return current seed.
        JMRNDM = HWRGET(SEED)

      ELSE

        WRITE(JMOUT,*) 'JMRNDM CALLED WITH ILLEGAL WOT=',WOT
        STOP

      ENDIF
      RETURN
      END

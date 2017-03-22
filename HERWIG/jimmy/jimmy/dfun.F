      DOUBLE PRECISION FUNCTION DFUN( N, xy )
C ----------------------------------------------------------------------------
C           Top level function call: integrated by DIVON4.
C           What function is actually called depends upon the value of
C           FN_TYPE in the JMFLAG common block.
C
C           3 = Eikonalising gamma-p xsec
C           4 = Calculating xsec for there being N (& only N) scatters.
C         101 = Calculating the total (uneikonalised) xsec.
C ----------------------------------------------------------------------------

      IMPLICIT NONE

#include "jimmy.inc"


      INTEGER N,i
      DOUBLE PRECISION XY(N)
      DOUBLE PRECISION JMKERN, EIKFUN, JMSNFN

      IF (FN_TYPE.EQ.101) THEN

C --    Uneikonalised cross section
        DFUN = JMKERN(N,XY)

      ELSE IF (FN_TYPE.EQ.3) THEN

C --    Eikonalising total xsec
        DFUN = EIKFUN(N,XY)

      ELSE IF (FN_TYPE.EQ.4) THEN

C --    Calculating xsec for N scatters
        DFUN = JMSNFN(N,XY)

      ELSE

C --    Illegal call.
        WRITE (JMOUT,*) 'DFUN:Illegal value of FN_TYPE!'
        STOP

      ENDIF

      RETURN
      END

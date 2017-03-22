*-- Author : Jon Butterworth
      DOUBLE PRECISION FUNCTION JMGAMZ( Z, WHICH, J )
C -------------------------------------------------------------------
C  Purpose: Returns total gamma-p cross section at given z, by
C           interpolating logarithmically between the points stored
C           during initialisation.
C           If the incoming particles are at fixed CM energy,
C           (e.g. pp) then input should be Z=-1.0, and J is set to 1
C           WHICH = 1 returns uneikonalised xsec
C           WHICH = 2 returns eikonalised xsec
C -------------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"


      DOUBLE PRECISION Z
      DOUBLE PRECISION ZSIN, LOG1, LOG2, Z1, Z2, LOGSIG, M, C
      INTEGER I, J, WHICH, LOC

      JMGAMZ = 0.D0

      IF (IERROR.NE.0) RETURN

      IF (WHICH.NE.1.AND.WHICH.NE.2) THEN
         WRITE(JMOUT,*) 'JMGAMZ called with illegal WHICH flag',WHICH
         STOP
      ENDIF

      LOC = WHICH+1

      IF (Z.LT.-0.5D0) THEN
        J=1
        JMGAMZ=JMARRY(LOC,J)
        RETURN
      ENDIF

      ZSIN = Z
      DO 1, I = 1, NPSIMP

        IF (ZSIN.LE.JMARRY( 1, I )) GOTO 2

 1    CONTINUE

      WRITE(JMOUT,*) 'JMGAMZ: Illegal Z input!', ZSIN, JMARRY(1,NPSIMP)
      IERROR = 200
      RETURN

C     We are interpolating between JMARRY(LOC,I-1) and JMARRY(LOC,I)
 2    CONTINUE

      J=I

C     Do a logarithmic interpolation if possible/sensible.
      IF (JMARRY(LOC,I-1).GT.0.AND.JMARRY(LOC,I).GT.0) THEN

        LOG1 = DLOG( JMARRY(LOC,I-1) )
        LOG2 = DLOG( JMARRY(LOC,I) )
        Z1 = JMARRY(1,I-1)
        Z2 = JMARRY(1,I)

        M = (LOG2-LOG1)/(Z2-Z1)
        C = LOG1-M*Z1

        LOGSIG = M*ZSIN + C

        JMGAMZ = EXP( LOGSIG )

      ELSE

C       Otherwise do a linear one.
        LOG1 = JMARRY(LOC,I-1)
        LOG2 = JMARRY(LOC,I)
        Z1 = JMARRY(1,I-1)
        Z2 = JMARRY(1,I)

        M = (LOG2-LOG1)/(Z2-Z1)
        C = LOG1-M*Z1

        LOGSIG = M*ZSIN + C

        JMGAMZ = LOGSIG

      ENDIF

      RETURN
      END

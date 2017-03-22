*CMZ :          17/07/95  11.08.55  by  Jonathan Butterworth
*-- Author : J. M. Butterworth
      SUBROUTINE JMPN( PROBS, ZDUB, K )
C --------------------------------------------------------------------------
C  Purpose: Returns the probabilities of N scatters at a given z, by
C           interpolating between the points stored during initialisation.
C           Returned arg. J is the lower edge of the s bin we are in.
C --------------------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"


      DOUBLE PRECISION ZDUB
      DOUBLE PRECISION Z, PROBS( MAXMS )
      DOUBLE PRECISION JMSC1, JMSC2, Z1, Z2, SCSIG, M, C
      INTEGER I, J, K

      Z=ZDUB
      DO I=1,MAXMS
         PROBS(I)=0.D0
      ENDDO

      IF (IERROR.NE.0) RETURN

      IF (ZDUB.LT.-0.5D0) THEN
        DO J=1,MAXMS
           PROBS(J)=JMARRY(4+J,1)
        ENDDO
        RETURN
      ENDIF

      DO 1, I=1, NPSIMP

        IF (Z.LE.JMARRY( 1, I )) GOTO 2

 1    CONTINUE

      WRITE(JMOUT,*) 'JMPN: ILLEGAL Z INPUT! ',Z
      IERROR=200
      RETURN

C     We are interpolating between JMARRY(J,I-1) and JMARRY(J,I)
C     Linear interpolation.
 2    CONTINUE

      K=I

      DO J=1, MAXMS

        IF (JMARRY(4+J,I-1).GT.0.0) THEN
          JMSC1 = JMARRY(4+J,I-1)
        ENDIF

        IF (JMARRY(4+J,I).GT.0.D0) THEN
          JMSC2 = JMARRY(4+J,I)
        ENDIF

        Z1 = JMARRY(1,I-1)
        Z2 = JMARRY(1,I)

        M = (JMSC2-JMSC1)/(Z2-Z1)
        C = JMSC1-M*Z1

        SCSIG = M*Z + C

        PROBS(J) = SCSIG

      ENDDO

      RETURN
      END

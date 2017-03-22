      SUBROUTINE JMKNIF(X1,X2,PT2,VETO)
C --------------------------------------------------------------------------
C Enforce virtual mass cutoffs of the structure function values
C generated.
C --------------------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"


      INTEGER I
      DOUBLE PRECISION HWBVMC
      DOUBLE PRECISION JMS, X1, X2, PT2, QL1, QL2
      LOGICAL VETO

      VETO = .FALSE.

      JMS=2.D0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2)

C     Demand that there's enough energy to produce at least a
C     pair of quarks. We can do this seperately for each subprocess
C     based upon their final state masses, but this will do for now.
      IF (PT2+RMASS(1)**2.GT.(X1*X2*YGAMMA*JMS)/4.D0) THEN
        VETO = .TRUE.
	IF (JMBUG.GT.3) THEN
          WRITE(JMOUT,8900) 'JMVETO: Not enough energy. Vetoing.'
	  WRITE(JMOUT,*) PT2,RMASS(1),(X1*X2*YGAMMA*JMS)/4.D0
	ENDIF
        RETURN
      ENDIF

C     Now apply the initial state cutoffs, process by process.
      QL1 = (1.D0-X1)*EMSCA
      QL2 = (1.D0-X2)*EMSCA

      DO I = 1, 13
        JMVETO(1,I) = 1.
        JMVETO(2,I) = 1.
        IF (QL1.LT.HWBVMC(I)) JMVETO(1,I)=0.
        IF (QL2.LT.HWBVMC(I)) JMVETO(2,I)=0.
      ENDDO

 8900 FORMAT(A)
      RETURN
      END




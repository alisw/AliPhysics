*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE COSLEG.FOR
*COPY COSLEG
*
*=== cosleg ===========================================================*
*
      FUNCTION COSLEG ( ACOLEG )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
      DIMENSION ACOLEG (0:2)
      REAL RNDM(2)
*
      IF ( ACOLEG (0) .NE. 0.5D+00 ) THEN
         ANORM = 0.5D+00 / ACOLEG (0)
         DO 100 I = 0,2
            ACOLEG (I) = ACOLEG (I) * ANORM
  100    CONTINUE
      END IF
      ACOF = 1.5D+00 * ACOLEG (2)
      BCOF = ACOLEG (1)
      CCOF = ACOLEG (0) - 0.5D+00 * ACOLEG (2)
      FORWRA = 0.3333333333333333D+00 * ACOF
      FORWRB = 0.5D+00 * BCOF
      FORWRC = CCOF
      FORWRD = FORWRA + FORWRB + FORWRC
      BACKWR = 1.D+00 - FORWRD
      CALL GRNDM(RNDM,1)
      RNDCOS = RNDM (1)
      IF ( RNDCOS .LT. FORWRD ) THEN
         RNDCOS = RNDCOS
         IF ( RNDCOS .LT. FORWRA ) THEN
            COSLEG = ( RNDCOS / FORWRA )**0.3333333333333333D+00
         ELSE IF ( RNDCOS .LT. FORWRA + FORWRB ) THEN
            COSLEG = SQRT ( ( RNDCOS - FORWRA ) / FORWRB )
         ELSE
            COSLEG = ( RNDCOS - FORWRA - FORWRB ) / FORWRC
         END IF
      ELSE
         COSLEG = - ( RNDCOS - FORWRD ) / BACKWR
         IF ( ACOF .LT. - ANGLGB ) THEN
            COSMAX = - 0.5D+00 * BCOF / ACOF
            IF ( COSMAX .GT. 0.D+00 ) THEN
               PROMAX = CCOF
            ELSE IF ( COSMAX .LT. -1.D+00 ) THEN
               PROMAX = ACOF - BCOF + CCOF
            ELSE
               PROMAX = ACOF * COSMAX**2 + BCOF * COSMAX + CCOF
            END IF
         ELSE IF ( ACOF .GT. ANGLGB ) THEN
            IF ( ACOF .GT. BCOF ) THEN
               PROMAX = ACOF - BCOF + CCOF
            ELSE
               PROMAX = CCOF
            END IF
         ELSE
            PROMAX = CCOF
         END IF
  200    CONTINUE
            FREJE = ( ACOF * COSLEG**2 + BCOF * COSLEG + CCOF ) / PROMAX
            CALL GRNDM(RNDM,2)
            IF ( RNDM (1) .GE. FREJE ) THEN
               COSLEG = - RNDM (2)
               GO TO 200
            END IF
      END IF
*=== End of subroutine cosleg =========================================*
      RETURN
      END

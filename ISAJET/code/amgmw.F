#include "isajet/pilot.h"
      FUNCTION AMGMW(I,J)
C
C          Get masses and widths from ISAJET commons for MadGraph
C          I = particle IDENT
C          J = 1 for mass
C            = 2 for width
C            = 3 for sin^2(theta)
C          Needed to avoid common block name clashes with MadGraph
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/wcon.inc"
#include "isajet/hcon.inc"
#include "isajet/sstype.inc"
      INTEGER I,J
      REAL AMGMW,AMASS
C
      IF(J.EQ.1) THEN
        AMGMW=AMASS(I)
      ELSEIF(J.EQ.2.AND.I.EQ.IDW) THEN
        AMGMW=WGAM(2)
      ELSEIF(J.EQ.2.AND.I.EQ.IDZ) THEN
        AMGMW=WGAM(4)
      ELSEIF(J.EQ.2.AND.I.EQ.IDH) THEN
        AMGMW=HGAM
      ELSEIF(J.EQ.3.AND.I.EQ.1) THEN
        AMGMW=SIN2W
      ELSE
        WRITE(ITLIS,*) 'ERROR IN AMGMW: I,J =',I,J
        STOP99
      ENDIF
      RETURN
      END

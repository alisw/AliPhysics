      DOUBLE PRECISION FUNCTION JMAREA( B2 )
c ----------------------------------------------------------------------
c   Purpose: Returns area overlap A(B)
c      Mods: 20-Aug-1993 JMB. Allow use of CERNLIB.
c      Mods: 27-Jan-1995  RW. Allow e+e- (not NAGLIB!)
c ----------------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"

      INTEGER NPHOTONS
      DOUBLE PRECISION B2,AREA
      DOUBLE PRECISION JMDBESK0, JMDBESK1, BESSEL(0:3)
      DOUBLE PRECISION X1, X2

      IF (B2.LE.0.D0) THEN
        IF (JMBUG.GT.11) WRITE(JMOUT,*)'**** JMAREA:B2=',B2
        JMAREA = 0.D0
        RETURN
      ENDIF
      X1=(JMV2*B2)**.5
      X2=(JMU2*B2)**.5

      NPHOTONS=0
      IF ((ABS(IDPDG(IPART1)).EQ.11).OR.(ABS(IDPDG(IPART1)).EQ.22))
     $     NPHOTONS=NPHOTONS+1
      IF ((ABS(IDPDG(IPART2)).EQ.11).OR.(ABS(IDPDG(IPART2)).EQ.22))
     $     NPHOTONS=NPHOTONS+1

      IF (NPHOTONS.EQ.1) THEN
C     `ep' type run
C     Modified Bessel function K0
         AREA = JMDBESK0(X1)
         AREA = (AREA-JMDBESK0(X2))*JMU2/(JMU2-JMV2)

C       Modified Bessel function K1
         AREA=(AREA-X2/TWO*JMDBESK1(X2))*JMU2*JMV2/(JMU2-JMV2)/
     &        (TWO*PIFAC)
         
      ELSE IF (NPHOTONS.EQ.2) THEN
         
C     `e+e-' type run.
C     Modified Bessel function K1
         AREA = JMV2*X1*JMDBESK1(X1)/(4.D0*PIFAC)
         
         
      ELSE IF (NPHOTONS.EQ.0) THEN
         
C     `pp' type run.
C     Modified Bessel function K3
         CALL JMDBSKA(X2,0,1,3,BESSEL)
         AREA = JMU2/96/PIFAC*(X2**3)*BESSEL(3)
         
      ENDIF

      JMAREA = AREA

      RETURN
      END

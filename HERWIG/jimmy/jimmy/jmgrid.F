      SUBROUTINE JMGRID
C     ------------------------------------------------------------------------
C     Purpose: Calculate the total hadronic cross section grid with
C     eikonalisation. This routine uses the functions QCDRES (to
C     calculate the uneikonalised gp cross section) and EIKON (to
C     eikonalise it).
C     ------------------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"

      INTEGER I
      DOUBLE PRECISION XSECN, EPFUN, H, Y, SUM, TMP, DELTA

      DELTA = 0.005D0

      IF (JCMVAR.GT.0) THEN

c     Step-size
         H = (ONE - JMZMIN)/FLOAT(NPSIMP)
         
c     Lower limit of integration
         Y = JMZMIN+DELTA
         
         SUM = EPFUN( Y, 0 )
         
         Y = JMZMIN
         
         DO 2, I=1, NPSIMP/2
            
            Y = Y + H
            TMP = EPFUN(Y, 2*I-1 )*FOUR
            SUM = SUM + TMP
            
            Y = Y + H
            TMP = EPFUN(Y, 2*I )*TWO
            SUM = SUM + TMP
            
 2       CONTINUE
         
C     Get the weight for the end point correct
         SUM = SUM - TMP/TWO
         XSECN = SUM*H/THREE

      ELSE
         
         SUM   = EPFUN( 1.D0, 1 )
         XSECN = SUM
         
      ENDIF
         
      RETURN
   50 FORMAT(8F10.5)
      END




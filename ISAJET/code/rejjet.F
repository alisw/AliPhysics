#include "isajet/pilot.h"
      LOGICAL FUNCTION REJJET() 
C-----------------------------------------------------------------------
C-
C-         This is called after EVOLVE for TWOJET and DRELLYAN events 
C-         to test the partons (jets). REJJET=.FALSE. keeps the event.
C-
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/jetset.inc"
      INTEGER I,IFLAV  
      REJJET=.FALSE.
C***************************************
C           Sample REJJET function which keeps the event if one of the 
C           outgoing partons is a heavy quark. 
C      REJJET=.TRUE. 
C      DO 1 I=1,NJSET    
C        IF(JDCAY(I).NE.0) GO TO 1   
C        IFLAV=IABS(JTYPE(I))    
C        IF(IFLAV.GE.4.AND.IFLAV.LT.9) THEN
C          REJJET=.FALSE.  
C          RETURN  
C        ENDIF
C  1   CONTINUE  
C***************************************
      RETURN    
      END   

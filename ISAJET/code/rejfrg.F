#include "isajet/pilot.h"
      LOGICAL FUNCTION REJFRG() 
C-----------------------------------------------------------------------
C-
C-         This is called after FRGMNT for TWOJET and DRELLYAN events 
C-         to test the fragmentation. REJFRG=.FALSE. keeps the event.
C-
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/partcl.inc"
      INTEGER I,ID 
      REAL PTL,CUTLEP,CUTNU 
      DATA CUTLEP/50./,CUTNU/100./  
      REJFRG=.FALSE.
C***************************************
C           Sample REJFRG function which keeps the event if it contains
C           any lepton satisfying
C             PT > CUTLEP (charged lepton)  
C             PT > CUTNU  (neutrino)    
C           Appropriate values of the cuts must be set by the user.
C      REJFRG=.TRUE. 
C      DO 1 I=1,NPTCL    
C        IF(IDCAY(I).NE.0) GO TO 1
C        ID=IABS(IDENT(I))   
C        IF(ID.LE.10.OR.ID.GE.20) GO TO 1    
C        PTL=SQRT(PPTCL(1,I)**2+PPTCL(2,I)**2)   
C        IF((ID.EQ.11.OR.ID.EQ.13.OR.ID.EQ.15).AND.PTL.GT.CUTNU) THEN
C          REJFRG=.FALSE.    
C          RETURN
C        ELSEIF((ID.EQ.12..OR.ID.EQ.14).AND.PTL.GT.CUTLEP) THEN 
C          REJFRG=.FALSE.    
C          RETURN
C        ENDIF
C   1  CONTINUE  
C***************************************
      RETURN    
      END   

* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYKCUT_HIJING(MCUT)   
    
C...Dummy routine, which the user can replace in order to make cuts on  
C...the kinematics on the parton level before the matrix elements are   
C...evaluated and the event is generated. The cross-section estimates   
C...will automatically take these cuts into account, so the given   
C...values are for the allowed phase space region only. MCUT=0 means    
C...that the event has passed the cuts, MCUT=1 that it has failed.  
#include "pypars_hijing.inc"
    
      MCUT=0    
    
      RETURN    
      END   

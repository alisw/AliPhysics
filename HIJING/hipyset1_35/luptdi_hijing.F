* $Id$
    
C*********************************************************************  
    
      SUBROUTINE LUPTDI_HIJING(KFL,PX,PY)  
    
C...Purpose: to generate transverse momentum according to a Gaussian.   
#include "ludat1_hijing.inc"
    
C...Generate p_T and azimuthal angle, gives p_x and p_y.    
      KFLA=IABS(KFL)    
      PT=PARJ(21)*SQRT(-LOG(MAX(1E-10,RLU_HIJING(0)))) 
      IF(MSTJ(91).EQ.1) PT=PARJ(22)*PT  
      IF(KFLA.EQ.0.AND.MSTJ(13).LE.0) PT=0. 
      PHI=PARU(2)*RLU_HIJING(0)    
      PX=PT*COS(PHI)    
      PY=PT*SIN(PHI)    
    
      RETURN    
      END   

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION RPLNML(RX,N,RA,MODE)                                             
      REAL RPLNML,RX,RA,RS                                                      
                                                                                
      DIMENSION RA(0:*)
                                                                                
      RS=0                                                                      
      IF(MODE .EQ. 1) THEN                                                      
       IF(N .GE. 0) RS=RA(N)                                                    
       DO 1 J = N-1,0,-1                                                        
       RS=RA(J)+RX*RS                                                           
    1  CONTINUE                                                                 
      ELSEIF(MODE .EQ. -1) THEN                                                 
       IF(N .GE. 0) RS=RA(0)                                                    
       DO 2 J = 1,N                                                             
       RS=RA(J)+RX*RS                                                           
    2  CONTINUE                                                                 
      ENDIF                                                                     
      RPLNML=RS                                                                 
      RETURN                                                                    
      END
                                                                                
      FUNCTION CPLNML(CX,N,CA,MODE)
      COMPLEX CPLNML,CX,CA,CS                                                   
                                                                                
      DIMENSION CA(0:*)                                                 

      CS=0                                                                      
      IF(MODE .EQ. 1) THEN                                                      
       IF(N .GE. 0) CS=CA(N)                                                    
       DO 3 J = N-1,0,-1                                                        
       CS=CA(J)+CX*CS                                                           
    3  CONTINUE                                                                 
      ELSEIF(MODE .EQ. -1) THEN                                                 
       IF(N .GE. 0) CS=CA(0)                                                    
       DO 4 J = 1,N                                                             
       CS=CA(J)+CX*CS                                                           
    4  CONTINUE                                                                 
      ENDIF                                                                     
      CPLNML=CS                                                                 
      RETURN                                                                    
      END                                                                       

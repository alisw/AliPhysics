*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if !defined(CERNLIB_DOUBLE)
      FUNCTION DPLNML(DX,N,DA,MODE)                                             
#include "gen/imp128.inc"
      CHARACTER*6 NAME
      NAME = 'DPLNML'
      CALL MTLPRT(NAME,'B105',
     +'not available on this machine - see documentation')
      DPLNML=0.
      RETURN
      END
#endif
#if defined(CERNLIB_DOUBLE)
      FUNCTION DPLNML(DX,N,DA,MODE)                                             
      DOUBLE PRECISION DPLNML,DX,DA,DS                                          
                                                                                
      DIMENSION DA(0:*)
                                                                                
      DS=0                                                                      
      IF(MODE .EQ. 1) THEN                                                      
       IF(N .GE. 0) DS=DA(N)                                                    
       DO 1 J = N-1,0,-1                                                        
       DS=DA(J)+DX*DS                                                           
    1  CONTINUE                                                                 
      ELSEIF(MODE .EQ. -1) THEN                                                 
       IF(N .GE. 0) DS=DA(0)                                                    
       DO 2 J = 1,N                                                             
       DS=DA(J)+DX*DS                                                           
    2  CONTINUE                                                                 
      ENDIF                                                                     
      DPLNML=DS                                                                 
      RETURN                                                                    
      END
                                                                                
      FUNCTION WPLNML(WX,N,WA,MODE)
      COMPLEX*16 WPLNML,WX,WA,WS                                                
                                                                                
      DIMENSION WA(0:*)                                                 

      WS=0                                                                      
      IF(MODE .EQ. 1) THEN                                                      
       IF(N .GE. 0) WS=WA(N)                                                    
       DO 3 J = N-1,0,-1                                                        
       WS=WA(J)+WX*WS                                                           
    3  CONTINUE                                                                 
      ELSEIF(MODE .EQ. -1) THEN                                                 
       IF(N .GE. 0) WS=WA(0)                                                    
       DO 4 J = 1,N                                                             
       WS=WA(J)+WX*WS                                                           
    4  CONTINUE                                                                 
      ENDIF                                                                     
      WPLNML=WS                                                                 
      RETURN                                                                    
      END                                                                       
#endif

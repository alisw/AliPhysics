*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION CPSIPG(Z,K)                                                      
      COMPLEX CPSIPG,Z                                                          
      COMPLEX*16 WPSIPG,W                                                       
      DOUBLE PRECISION D                                                        
                                                                                
      SROUND(D)=D+(D-SNGL(D))                                                   
      W=Z                                                                       
      W=WPSIPG(W,K)                                                             
      CPSIPG=CMPLX(SROUND(DREAL(W)),SROUND(DIMAG(W)))                           
      RETURN                                                                    
      END                                                                       

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION WPSIPG(Z,K)                                                      
                                                                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      COMPLEX*16 WPSIPG,Z,U,V,H,R,P,GCMPLX                                      
      CHARACTER NAME*(*)                                                        
      CHARACTER*80 ERRTXT                                                       
      PARAMETER (NAME = 'CPSIPG/WPSIPG')                                        
      DIMENSION C(6,0:4),FCT(-1:4),SGN(0:4)                                     
                                                                                
      PARAMETER (DELTA = 5D-13)                                                 
      PARAMETER (R1 = 1, HF = R1/2)                                             
      PARAMETER (PI = 3.14159 26535 89793 24D0)                                 
      PARAMETER (C1 = PI**2, C2 = 2*PI**3, C3 = 2*PI**4, C4 = 8*PI**5)          
                                                                                
      DATA SGN /-1,+1,-1,+1,-1/, FCT /0,1,1,2,6,24/                             
                                                                                
      DATA C(1,0) / 8.33333 33333 33333 33D-2/                                  
      DATA C(2,0) /-8.33333 33333 33333 33D-3/                                  
      DATA C(3,0) / 3.96825 39682 53968 25D-3/                                  
      DATA C(4,0) /-4.16666 66666 66666 67D-3/                                  
      DATA C(5,0) / 7.57575 75757 57575 76D-3/                                  
      DATA C(6,0) /-2.10927 96092 79609 28D-2/                                  
                                                                                
      DATA C(1,1) / 1.66666 66666 66666 67D-1/                                  
      DATA C(2,1) /-3.33333 33333 33333 33D-2/                                  
      DATA C(3,1) / 2.38095 23809 52380 95D-2/                                  
      DATA C(4,1) /-3.33333 33333 33333 33D-2/                                  
      DATA C(5,1) / 7.57575 75757 57575 76D-2/                                  
      DATA C(6,1) /-2.53113 55311 35531 14D-1/                                  
                                                                                
      DATA C(1,2) / 5.00000 00000 00000 00D-1/                                  
      DATA C(2,2) /-1.66666 66666 66666 67D-1/                                  
      DATA C(3,2) / 1.66666 66666 66666 67D-1/                                  
      DATA C(4,2) /-3.00000 00000 00000 00D-1/                                  
      DATA C(5,2) / 8.33333 33333 33333 33D-1/                                  
      DATA C(6,2) /-3.29047 61904 76190 48D+0/                                  
                                                                                
      DATA C(1,3) / 2.00000 00000 00000 00D+0/                                  
      DATA C(2,3) /-1.00000 00000 00000 00D+0/                                  
      DATA C(3,3) / 1.33333 33333 33333 33D+0/                                  
      DATA C(4,3) /-3.00000 00000 00000 00D+0/                                  
      DATA C(5,3) / 1.00000 00000 00000 00D+1/                                  
      DATA C(6,3) /-4.60666 66666 66666 67D+1/                                  
                                                                                
      DATA (C(I,4),I=1,6) / 10, -7, 12, -33, 130, -691/                         
                                                                                
      GCMPLX(X,Y)=DCMPLX(X,Y)                                                   
                                                                                
      U=Z                                                                       
      X=U                                                                       
      A=ABS(X)                                                                  
      IF(K .LT. 0 .OR. K .GT. 4) THEN                                           
       H=0                                                                      
       WRITE(ERRTXT,101) K                                                      
       CALL MTLPRT(NAME,'C317.1',ERRTXT)                                        
      ELSEIF(ABS(IMAG(U)) .LT. DELTA .AND. ABS(X+NINT(A)) .LT. DELTA)           
     1                                                        THEN              
       H=0                                                                      
       WRITE(ERRTXT,102) X                                                      
       CALL MTLPRT(NAME,'C317.2',ERRTXT)                                        
      ELSE                                                                      
       K1=K+1                                                                   
       IF(X .LT. 0) U=-U                                                        
       V=U                                                                      
       H=0                                                                      
       IF(A .LT. 15) THEN                                                       
        H=1/V**K1                                                               
        DO 1 I = 1,14-INT(A)                                                    
        V=V+1                                                                   
    1   H=H+1/V**K1                                                             
        V=V+1                                                                   
       END IF                                                                   
       R=1/V**2                                                                 
       P=R*C(6,K)                                                               
       DO 2 I = 5,1,-1                                                          
    2  P=R*(C(I,K)+P)                                                           
       H=SGN(K)*(FCT(K)*H+(V*(FCT(K-1)+P)+HF*FCT(K))/V**K1)                     
       IF(K .EQ. 0) H=H+LOG(V)                                                  
       IF(X .LT. 0) THEN                                                        
        V=PI*U                                                                  
        X=V                                                                     
        Y=IMAG(V)                                                               
        A=SIN(X)                                                                
        B=COS(X)                                                                
        T=TANH(Y)                                                               
        P=GCMPLX(B,-A*T)/GCMPLX(A,B*T)                                          
        IF(K .EQ. 0) THEN                                                       
         H=H+1/U+PI*P                                                           
        ELSEIF(K .EQ. 1) THEN                                                   
         H=-H+1/U**2+C1*(P**2+1)                                                
        ELSEIF(K .EQ. 2) THEN                                                   
         H=H+2/U**3+C2*P*(P**2+1)                                               
        ELSEIF(K .EQ. 3) THEN                                                   
         R=P**2                                                                 
         H=-H+6/U**4+C3*((3*R+4)*R+1)                                           
        ELSEIF(K .EQ. 4) THEN                                                   
         R=P**2                                                                 
         H=H+24/U**5+C4*P*((3*R+5)*R+2)                                         
        ENDIF                                                                   
       ENDIF                                                                    
      ENDIF                                                                     
      WPSIPG=H                                                                  
      RETURN                                                                    
  101 FORMAT('K = ',I5,'  (< 0  OR  > 4)')                                      
  102 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',F8.1)                    
      END                                                                       

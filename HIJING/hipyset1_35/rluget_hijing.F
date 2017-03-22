* $Id$
    
C*********************************************************************  
    
      SUBROUTINE RLUGET_HIJING(LFN,MOVE)   
    
C...Purpose: to dump the state of the random number generator on a file 
C...for subsequent startup from this state onwards. 
#include "ludatr_hijing.inc"
      CHARACTER CHERR*8 
    
C...Backspace required number of records (or as many as there are). 
      IF(MOVE.LT.0) THEN    
        NBCK=MIN(MRLU(6),-MOVE) 
        DO 100 IBCK=1,NBCK  
  100   BACKSPACE(LFN,ERR=110,IOSTAT=IERR)  
        MRLU(6)=MRLU(6)-NBCK    
      ENDIF 
    
C...Unformatted write on unit LFN.  
      WRITE(LFN,ERR=110,IOSTAT=IERR) (MRLU(I1),I1=1,5), 
     &(RRLU(I2),I2=1,100)   
      MRLU(6)=MRLU(6)+1 
      RETURN    
    
C...Write error.    
  110 WRITE(CHERR,'(I8)') IERR  
      CALL LUERRM_HIJING(18
     $     ,'(RLUGET_HIJING:) error when accessing file, IOSTAT ='/
     $     /CHERR)    
    
      RETURN    
      END   

* $Id$
    
C*********************************************************************  
    
      SUBROUTINE LUGIVE_HIJING(CHIN)   
    
C...Purpose: to set values of commonblock variables.    
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "ludat3_hijing.inc"
#include "ludat4_hijing.inc"
      CHARACTER CHIN*(*),CHFIX*104,CHBIT*104,CHOLD*8,CHNEW*8,   
     &CHNAM*4,CHVAR(17)*4,CHALP(2)*26,CHIND*8,CHINI*10,CHINR*16 
      DATA CHVAR/'N','K','P','V','MSTU','PARU','MSTJ','PARJ','KCHG',    
     &'PMAS','PARF','VCKM','MDCY','MDME','BRAT','KFDP','CHAF'/  
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',  
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/ 
    
C...Length of character variable. Subdivide it into instructions.   
      IF(MSTU(12).GE.1) CALL LULIST_HIJING(0)  
      CHBIT=CHIN//' '   
      LBIT=101  
  100 LBIT=LBIT-1   
      IF(CHBIT(LBIT:LBIT).EQ.' ') GOTO 100  
      LTOT=0    
      DO 110 LCOM=1,LBIT    
      IF(CHBIT(LCOM:LCOM).EQ.' ') GOTO 110  
      LTOT=LTOT+1   
      CHFIX(LTOT:LTOT)=CHBIT(LCOM:LCOM) 
  110 CONTINUE  
      LLOW=0    
  120 LHIG=LLOW+1   
  130 LHIG=LHIG+1   
      IF(LHIG.LE.LTOT.AND.CHFIX(LHIG:LHIG).NE.';') GOTO 130 
      LBIT=LHIG-LLOW-1  
      CHBIT(1:LBIT)=CHFIX(LLOW+1:LHIG-1)    
    
C...Identify commonblock variable.  
      LNAM=1    
  140 LNAM=LNAM+1   
      IF(CHBIT(LNAM:LNAM).NE.'('.AND.CHBIT(LNAM:LNAM).NE.'='.AND.   
     &LNAM.LE.4) GOTO 140   
      CHNAM=CHBIT(1:LNAM-1)//' '    
      DO 150 LCOM=1,LNAM-1  
      DO 150 LALP=1,26  
  150 IF(CHNAM(LCOM:LCOM).EQ.CHALP(1)(LALP:LALP)) CHNAM(LCOM:LCOM)= 
     &CHALP(2)(LALP:LALP)   
      IVAR=0    
      DO 160 IV=1,17    
  160 IF(CHNAM.EQ.CHVAR(IV)) IVAR=IV    
      IF(IVAR.EQ.0) THEN    
         CALL LUERRM_HIJING(18
     $        ,'(LUGIVE_HIJING:) do not recognize variable '//CHNAM)   
        LLOW=LHIG   
        IF(LLOW.LT.LTOT) GOTO 120   
        RETURN  
      ENDIF 
    
C...Identify any indices.   
      I=0   
      J=0   
      IF(CHBIT(LNAM:LNAM).EQ.'(') THEN  
        LIND=LNAM   
  170   LIND=LIND+1 
        IF(CHBIT(LIND:LIND).NE.')'.AND.CHBIT(LIND:LIND).NE.',') GOTO 170    
        CHIND=' '   
        IF((CHBIT(LNAM+1:LNAM+1).EQ.'C'.OR.CHBIT(LNAM+1:LNAM+1).EQ.'c').    
     &  AND.(IVAR.EQ.9.OR.IVAR.EQ.10.OR.IVAR.EQ.13.OR.IVAR.EQ.17)) THEN 
          CHIND(LNAM-LIND+11:8)=CHBIT(LNAM+2:LIND-1)    
          READ(CHIND,'(I8)') I1 
          I=LUCOMP_HIJING(I1)  
        ELSE    
          CHIND(LNAM-LIND+10:8)=CHBIT(LNAM+1:LIND-1)    
          READ(CHIND,'(I8)') I  
        ENDIF   
        LNAM=LIND   
        IF(CHBIT(LNAM:LNAM).EQ.')') LNAM=LNAM+1 
      ENDIF 
      IF(CHBIT(LNAM:LNAM).EQ.',') THEN  
        LIND=LNAM   
  180   LIND=LIND+1 
        IF(CHBIT(LIND:LIND).NE.')'.AND.CHBIT(LIND:LIND).NE.',') GOTO 180    
        CHIND=' '   
        CHIND(LNAM-LIND+10:8)=CHBIT(LNAM+1:LIND-1)  
        READ(CHIND,'(I8)') J    
        LNAM=LIND+1 
      ENDIF 
    
C...Check that indices allowed and save old value.  
      IERR=1    
      IF(CHBIT(LNAM:LNAM).NE.'=') GOTO 190  
      IF(IVAR.EQ.1) THEN    
        IF(I.NE.0.OR.J.NE.0) GOTO 190   
        IOLD=N  
      ELSEIF(IVAR.EQ.2) THEN    
        IF(I.LT.1.OR.I.GT.MSTU(4).OR.J.LT.1.OR.J.GT.5) GOTO 190 
        IOLD=K(I,J) 
      ELSEIF(IVAR.EQ.3) THEN    
        IF(I.LT.1.OR.I.GT.MSTU(4).OR.J.LT.1.OR.J.GT.5) GOTO 190 
        ROLD=P(I,J) 
      ELSEIF(IVAR.EQ.4) THEN    
        IF(I.LT.1.OR.I.GT.MSTU(4).OR.J.LT.1.OR.J.GT.5) GOTO 190 
        ROLD=V(I,J) 
      ELSEIF(IVAR.EQ.5) THEN    
        IF(I.LT.1.OR.I.GT.200.OR.J.NE.0) GOTO 190   
        IOLD=MSTU(I)    
      ELSEIF(IVAR.EQ.6) THEN    
        IF(I.LT.1.OR.I.GT.200.OR.J.NE.0) GOTO 190   
        ROLD=PARU(I)    
      ELSEIF(IVAR.EQ.7) THEN    
        IF(I.LT.1.OR.I.GT.200.OR.J.NE.0) GOTO 190   
        IOLD=MSTJ(I)    
      ELSEIF(IVAR.EQ.8) THEN    
        IF(I.LT.1.OR.I.GT.200.OR.J.NE.0) GOTO 190   
        ROLD=PARJ(I)    
      ELSEIF(IVAR.EQ.9) THEN    
        IF(I.LT.1.OR.I.GT.MSTU(6).OR.J.LT.1.OR.J.GT.3) GOTO 190 
        IOLD=KCHG(I,J)  
      ELSEIF(IVAR.EQ.10) THEN   
        IF(I.LT.1.OR.I.GT.MSTU(6).OR.J.LT.1.OR.J.GT.4) GOTO 190 
        ROLD=PMAS(I,J)  
      ELSEIF(IVAR.EQ.11) THEN   
        IF(I.LT.1.OR.I.GT.2000.OR.J.NE.0) GOTO 190  
        ROLD=PARF(I)    
      ELSEIF(IVAR.EQ.12) THEN   
        IF(I.LT.1.OR.I.GT.4.OR.J.LT.1.OR.J.GT.4) GOTO 190   
        ROLD=VCKM(I,J)  
      ELSEIF(IVAR.EQ.13) THEN   
        IF(I.LT.1.OR.I.GT.MSTU(6).OR.J.LT.1.OR.J.GT.3) GOTO 190 
        IOLD=MDCY(I,J)  
      ELSEIF(IVAR.EQ.14) THEN   
        IF(I.LT.1.OR.I.GT.MSTU(7).OR.J.LT.1.OR.J.GT.2) GOTO 190 
        IOLD=MDME(I,J)  
      ELSEIF(IVAR.EQ.15) THEN   
        IF(I.LT.1.OR.I.GT.MSTU(7).OR.J.NE.0) GOTO 190   
        ROLD=BRAT(I)    
      ELSEIF(IVAR.EQ.16) THEN   
        IF(I.LT.1.OR.I.GT.MSTU(7).OR.J.LT.1.OR.J.GT.5) GOTO 190 
        IOLD=KFDP(I,J)  
      ELSEIF(IVAR.EQ.17) THEN   
        IF(I.LT.1.OR.I.GT.MSTU(6).OR.J.NE.0) GOTO 190   
        CHOLD=CHAF(I)   
      ENDIF 
      IERR=0    
  190 IF(IERR.EQ.1) THEN    
         CALL LUERRM_HIJING(18,'(LUGIVE_HIJING:) unallowed indices for '
     $        //CHBIT(1:LNAM-1))    
        LLOW=LHIG   
        IF(LLOW.LT.LTOT) GOTO 120   
        RETURN  
      ENDIF 
    
C...Print current value of variable. Loop back. 
      IF(LNAM.GE.LBIT) THEN 
        CHBIT(LNAM:14)=' '  
        CHBIT(15:60)=' has the value                                '   
        IF(IVAR.EQ.1.OR.IVAR.EQ.2.OR.IVAR.EQ.5.OR.IVAR.EQ.7.OR. 
     &  IVAR.EQ.9.OR.IVAR.EQ.13.OR.IVAR.EQ.14.OR.IVAR.EQ.16) THEN   
          WRITE(CHBIT(51:60),'(I10)') IOLD  
        ELSEIF(IVAR.NE.17) THEN 
          WRITE(CHBIT(47:60),'(F14.5)') ROLD    
        ELSE    
          CHBIT(53:60)=CHOLD    
        ENDIF   
        IF(MSTU(13).GE.1) WRITE(MSTU(11),1000) CHBIT(1:60)  
        LLOW=LHIG   
        IF(LLOW.LT.LTOT) GOTO 120   
        RETURN  
      ENDIF 
    
C...Read in new variable value. 
      IF(IVAR.EQ.1.OR.IVAR.EQ.2.OR.IVAR.EQ.5.OR.IVAR.EQ.7.OR.   
     &IVAR.EQ.9.OR.IVAR.EQ.13.OR.IVAR.EQ.14.OR.IVAR.EQ.16) THEN 
        CHINI=' '   
        CHINI(LNAM-LBIT+11:10)=CHBIT(LNAM+1:LBIT)   
        READ(CHINI,'(I10)') INEW    
      ELSEIF(IVAR.NE.17) THEN   
        CHINR=' '   
        CHINR(LNAM-LBIT+17:16)=CHBIT(LNAM+1:LBIT)   
        READ(CHINR,'(F16.2)') RNEW  
      ELSE  
        CHNEW=CHBIT(LNAM+1:LBIT)//' '   
      ENDIF 
    
C...Store new variable value.   
      IF(IVAR.EQ.1) THEN    
        N=INEW  
      ELSEIF(IVAR.EQ.2) THEN    
        K(I,J)=INEW 
      ELSEIF(IVAR.EQ.3) THEN    
        P(I,J)=RNEW 
      ELSEIF(IVAR.EQ.4) THEN    
        V(I,J)=RNEW 
      ELSEIF(IVAR.EQ.5) THEN    
        MSTU(I)=INEW    
      ELSEIF(IVAR.EQ.6) THEN    
        PARU(I)=RNEW    
      ELSEIF(IVAR.EQ.7) THEN    
        MSTJ(I)=INEW    
      ELSEIF(IVAR.EQ.8) THEN    
        PARJ(I)=RNEW    
      ELSEIF(IVAR.EQ.9) THEN    
        KCHG(I,J)=INEW  
      ELSEIF(IVAR.EQ.10) THEN   
        PMAS(I,J)=RNEW  
      ELSEIF(IVAR.EQ.11) THEN   
        PARF(I)=RNEW    
      ELSEIF(IVAR.EQ.12) THEN   
        VCKM(I,J)=RNEW  
      ELSEIF(IVAR.EQ.13) THEN   
        MDCY(I,J)=INEW  
      ELSEIF(IVAR.EQ.14) THEN   
        MDME(I,J)=INEW  
      ELSEIF(IVAR.EQ.15) THEN   
        BRAT(I)=RNEW    
      ELSEIF(IVAR.EQ.16) THEN   
        KFDP(I,J)=INEW  
      ELSEIF(IVAR.EQ.17) THEN   
        CHAF(I)=CHNEW   
      ENDIF 
    
C...Write old and new value. Loop back. 
      CHBIT(LNAM:14)=' '    
      CHBIT(15:60)=' changed from                to               ' 
      IF(IVAR.EQ.1.OR.IVAR.EQ.2.OR.IVAR.EQ.5.OR.IVAR.EQ.7.OR.   
     &IVAR.EQ.9.OR.IVAR.EQ.13.OR.IVAR.EQ.14.OR.IVAR.EQ.16) THEN 
        WRITE(CHBIT(33:42),'(I10)') IOLD    
        WRITE(CHBIT(51:60),'(I10)') INEW    
      ELSEIF(IVAR.NE.17) THEN   
        WRITE(CHBIT(29:42),'(F14.5)') ROLD  
        WRITE(CHBIT(47:60),'(F14.5)') RNEW  
      ELSE  
        CHBIT(35:42)=CHOLD  
        CHBIT(53:60)=CHNEW  
      ENDIF 
      IF(MSTU(13).GE.1) WRITE(MSTU(11),1000) CHBIT(1:60)    
      LLOW=LHIG 
      IF(LLOW.LT.LTOT) GOTO 120 
    
C...Format statement for output on unit MSTU(11) (by default 6).    
 1000 FORMAT(5X,A60)    
    
      RETURN    
      END   

#include "isajet/pilot.h"
      SUBROUTINE IRMOV0
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-      remove 0's from PJSET
C-
C-   Created  15-OCT-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
      INTEGER NCOUNT,I,J,K
C----------------------------------------------------------------------
C
C         remove zeroes
      NCOUNT=NJSET
      DO 160 I=3,NJSET  
  151   IF (PJSET(4,I).EQ.0.AND.I.LT.NCOUNT) THEN
          DO 155 K=I+1,NCOUNT
            DO 154 J=1,5    
              PJSET(J,K-1)=PJSET(J,K)   
  154       CONTINUE    
            JORIG(K-1)=JORIG(K) 
            JTYPE(K-1)=JTYPE(K) 
            JDCAY(K-1)=JDCAY(K) 
            ZZC(K-1)=ZZC(K)
            JMATCH(K-1)=JMATCH(K)
            IF(JMATCH(K-1).GT.I) JMATCH(K-1)=JMATCH(K-1)-1
  155     CONTINUE  
          NCOUNT=NCOUNT-1
          GOTO 151
        ENDIF
  160 CONTINUE  
      NJSET=NCOUNT  
C          remove last one if 0
      IF(PJSET(4,NJSET).EQ.0) NJSET=NJSET-1  
  999 RETURN
      END

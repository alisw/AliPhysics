*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:40  mclareni
* Kernlib
*
*
      SUBROUTINE JUMPX2 (IPARA,IPARB)
C
C CERN PROGLIB# Z042    JUMPX2          .VERSION KERNMPW  1.00  891208
C ORIG. 21/04/88 FCA
C
C-    To transfer to the user routine TARGET (say) with 2 parameters
C-    three steps are needed :
 
C- 1) EXTERNAL TARGET              to get the address of TARGET
C-    IADR = JUMPAD (TARGET)
 
C- 2) CALL JUMPST (IADR)           to set the tranfer address
 
C- 3) CALL JUMPX2 (par1,par2)      to transfer
 
      DIMENSION IPARA(9), IPARB(9)
 
      SAVE      JTARGE
 
      CALL JUMPY2 (%VAL(JTARGE),IPARA,IPARB)
      RETURN
 
      ENTRY JUMPX1 (IPARA)
 
      CALL JUMPY1 (%VAL(JTARGE),IPARA)
      RETURN
 
      ENTRY JUMPX0
 
      CALL JUMPY0 (%VAL(JTARGE))
      RETURN
 
      ENTRY JUMPST (JADR)
 
      JTARGE = JADR
      END
C-----------------------------------------------------
      SUBROUTINE JUMPY2 (TARSUB,IPARA,IPARB)
      DIMENSION IPARA(9), IPARB(9)
      EXTERNAL  TARSUB
 
      CALL TARSUB (IPARA,IPARB)
      END
C-----------------------------------------------------
      SUBROUTINE JUMPY1 (TARSUB,IPARA)
      DIMENSION IPARA(9)
      EXTERNAL  TARSUB
 
      CALL TARSUB (IPARA)
      END
C-----------------------------------------------------
      SUBROUTINE JUMPY0 (TARSUB)
      EXTERNAL  TARSUB
 
      CALL TARSUB
      END

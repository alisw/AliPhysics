*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:41  mclareni
* Kernlib
*
*
      SUBROUTINE JUMPX2 (IPARA,IPARB)
C
C CERN PROGLIB# Z042    JUMPX2          .VERSION KERNVAX  2.32  890926
C ORIG. 21/04/88 JZ
C MOD.  12/09/89 JZ, add entries JUMPX1, JUMPX0
C
C-    To transfer to the user routine TARGET (say) with 2 parameters
C-    three steps are needed :

C- 1) EXTERNAL TARGET              to get the address of TARGET
C-    IADR = JUMPAD (TARGET)

C- 2) CALL JUMPST (IADR)           to set the tranfer address

C- 3) CALL JUMPX2 (par1,par2)      to transfer

      DIMENSION IPARA(9), IPARB(9)
      SAVE JPOINT

      CALL  JUMPY2 (%VAL(JPOINT),IPARA,IPARB)
      RETURN

      ENTRY JUMPX1 (IPARA)
      CALL  JUMPY1 (%VAL(JPOINT),IPARA)
      RETURN

      ENTRY JUMPX0
      CALL  JUMPY0 (%VAL(JPOINT))
      RETURN

C-    to set the tranfer address
      ENTRY JUMPST (JADR)
      JPOINT = JADR
      RETURN
      END
C-----------------------------------------------------
      SUBROUTINE JUMPY2 (TARSUB,IPARA,IPARB)
      DIMENSION IPARA(9), IPARB(9)
      EXTERNAL  TARSUB

      CALL  TARSUB (IPARA,IPARB)
      RETURN

      ENTRY JUMPY1 (TARSUB,IPARA)
      CALL  TARSUB (IPARA)
      RETURN

      ENTRY JUMPY0 (TARSUB)
      CALL  TARSUB
      RETURN
      END

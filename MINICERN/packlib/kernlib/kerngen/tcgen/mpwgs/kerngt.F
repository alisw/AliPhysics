*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
!!S SEG_KERNGT
      SUBROUTINE KERNGT(LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNMPW  1.00  891208
C ORIG. 12/12/83  JZ
C
C--   STORE PAM TITLE INTO /SLATE/ FOR PRINTING BY GSHOWF
C
      COMMON /SLATE/ CH(40)
C
      ENTRY      GSHOW (LUN)
C
      CALL UBLOW (
     + 40H          KERNMPW  1.02  921009 00.35                 
     +,  CH,40)
      CALL GSHOWF (LUN)
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNSGI  1.07  940120
C ORIG. 16/03/89  JZ
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF
C
      COMMON /SLATE/ MCH(40)
      CHARACTER    VIDQQ*(*)
      PARAMETER   (VIDQQ = '@(#)' //
     + ' KERNSGI  1.05  940511 18.45           >')

      CALL CLEFT (VIDQQ,1,0)
      CALL UCTOH1 (
     + ' KERNSGI  1.05  940511 18.45           '
     +,  MCH,40)
      CALL GSHOWF (LUN)
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:57  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNIRT  1.08  940531
C ORIG. 12/12/83  JZ
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF
C
      COMMON /SLATE/ MCH(40)
      CHARACTER    VIDQQ*(*)
      PARAMETER   (VIDQQ = '@(#)' //
     + ' KERNIRT  1.08  940526 11.05            >')

      CALL CLEFT (VIDQQ,1,0)
      CALL UCTOH1 (
     + ' KERNIRT  1.08  940526 11.05            '
     +,  MCH,40)
      CALL GSHOWF (LUN)
      RETURN
      END

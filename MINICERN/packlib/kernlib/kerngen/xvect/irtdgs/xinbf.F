*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:36  mclareni
* Kernlib
*
*
      SUBROUTINE XINBF (LUN,XBV,NXB)
C
C CERN PROGLIB# Z203    XINBF           .VERSION KERNIRT  1.06  930811
C
      DIMENSION    XBV(NXB)
C
      READ (LUN,IOSTAT=JSTAT,NUM=NBYTES) XBV
      IF   (JSTAT)                 61, 24, 64
   24 NXB = MAX (1, NBYTES/8)
      RETURN
C
   61 NXB = 0
      RETURN
C
   64 NXB = -JSTAT
      READ (LUN) N
      RETURN
      END

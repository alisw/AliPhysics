*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:36  mclareni
* Kernlib
*
*
      SUBROUTINE XINBS (LUN,XAV,NA,XBV,NXB)
C
C CERN PROGLIB# Z203    XINBS           .VERSION KERNIRT  1.01  900319
C
      COMMON /SLATE/ NXR, DUMMY(39)
      DIMENSION    XBV(NXB), XAV(NA)
C
      READ (LUN,IOSTAT=JSTAT,NUM=NBYTES) NR,XAV,XBV
      IF   (JSTAT)                 61, 24, 64
   24 NXR = NR
      NR  = MIN (NR,NXB)
      NXB = MAX (NR,1)
      RETURN
C
   61 NXB = 0
      RETURN
C
   64 NXB = -JSTAT
      READ (LUN) N
      RETURN
      END

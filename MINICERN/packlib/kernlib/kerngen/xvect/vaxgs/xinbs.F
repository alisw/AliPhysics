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
C CERN PROGLIB# Z203    XINBS           .VERSION KERNVAX  2.29  880809
C ORIG.  8/02/88  JZ
C
      COMMON /SLATE/ NXR, DUMMY(39)
      DIMENSION    XBV(99), XAV(NA)

      READ (LUN,IOSTAT=JSTAT) NR,XAV,(XBV(J),J=1,MIN(NR,NXB))
      IF   (JSTAT)           61, 24, 64
   24 NXR = NR
      NR  = MIN (NR,NXB)
      NXB = MAX (NR,1)
      RETURN

   61 NXB = 0
      RETURN

   64 NXB = -JSTAT
      RETURN
      END

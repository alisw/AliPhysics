*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:36  mclareni
* Kernlib
*
*
      SUBROUTINE XINB  (LUN,XBV,NXB)
C
C CERN PROGLIB# Z203    XINB            .VERSION KERNVAX  2.29  880809
C ORIG.  8/02/88  JZ
C
      COMMON /SLATE/ NXR, DUMMY(39)
      DIMENSION    XBV(99)

      READ (LUN,IOSTAT=JSTAT) NR, (XBV(J),J=1,MIN(NR,NXB))
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

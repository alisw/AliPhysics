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
C CERN PROGLIB# Z203    XINB            .VERSION KERNMPW  1.00  891208
C ORIG. 06/02/89 JZ
C
      COMMON /SLATE/ NXR, DUMMY(39)
      INTEGER      XBV(999), MBUF(512)
 
      READ (LUN,END=61,ERR=64) NR,(XBV(J),J=1,MIN(NR,NXB))
      NXR = NR
      NR  = MIN(NR,NXB)
      NSK = NXR - NR
      NXB = MAX(1,NR)
 
   24 IF (NSK.LE.0)          RETURN
      NR = MIN (NSK,512)
      READ (LUN,END=61,ERR=64) (MBUF(J),J=1,NR)
      NSK = NSK - NR
      GO TO 24
 
   61 NXB = 0
      RETURN
 
   64 NXB = -7
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      SUBROUTINE PKBYT (MIV,MBV,JTHP,NINTP,NBITS)
C
C CERN PROGLIB# M422    PKBYT           .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C
      DIMENSION    MIV(99), MBV(99), JTHP(9), NINTP(9), NBITS(2)


      JTH  = JTHP(1)
      NINT = NINTP(1)
      IF (NINT.LE.0)         RETURN

      NZB  = NBITS(1)
      IF (NZB.GT.0)          GO TO 11
      NZB  = 1
      NPWD = 32
      MASK = 1
      GO TO 12

   11 NPWD = NBITS(2)
#if defined(CERNLIB_BUGLRSHFT)
      MASK = ishft (-1, NZB-32)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MASK = lrshft (-1, 32-NZB)
#endif

   12 JBV  = 1
      JIV  = 0
      IF (JTH.LT.2)          GO TO 21
      JBV  = (JTH-1)/NPWD + 1
      JPOS = JTH - (JBV-1)*NPWD - 1
      IF (JPOS.EQ.0)         GO TO 21
      NL   = JPOS*NZB
      MSKU = LSHIFT (MASK,NL)
      JIVE = NPWD - JPOS
      GO TO 22

C--                PACK EACH WORD

   21 NL   = 0
      MSKU = MASK
      JIVE = JIV + NPWD
   22 JIVE = MIN0 (NINT,JIVE)
      IZW  = MBV(JBV)

   24 JIV  = JIV + 1
      IZW  = OR (AND(NOT(MSKU),IZW),
     +           AND(MSKU,LSHIFT(MIV(JIV),NL)))
      IF (JIV.EQ.JIVE)       GO TO 27
      NL   = NL + NZB
      MSKU = LSHIFT (MSKU,NZB)
      GO TO 24

   27 MBV(JBV) = IZW
      IF (JIV.EQ.NINT)       RETURN
      JBV  = JBV + 1
      GO TO 21
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fziphr32.inc
*
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
      IF (IACMOI.EQ.0)             GO TO 36
      NWR = 4*NWMIN
      CALL IOREAD (LUNI,LQ(LIN),NWR,ISW)
      IF (ISW.EQ.0)                GO TO 39
      IF (ISW.EQ.1)                GO TO 843
      IF (ISW.EQ.2)                GO TO 841
      IF (ISW.EQ.3)                GO TO 843
      GO TO 841
#endif
*
* fziphr32.inc
*
#if defined(CERNLIB_QMND3)
      IF (IACMOI.EQ.0)             GO TO 36
      ISW = MAGTP (0,LQ(LIN),LUNI,NWMIN,NWR)
      IF (ISW.EQ.0)                GO TO 39
      IF (ISW.EQ.3)                GO TO 841
      GO TO 843
*    --------------  conversions for input   -------------------
#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzophr44.inc
*
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
      IF (IACMOX.EQ.0)             GO TO 48
      CALL IORITE (LUNX,LQ(LWR),4*NWR,ISW)
      IF (ISW.EQ.0)                GO TO 999
      IF (ISW.NE.2)                GO TO 45
      IF (IQ(KQSP+LQFX+37)+IQ(KQSP+LQFX+38).EQ.0)  GO TO 45
      IQ(KQSP+LQFX+37) = IQ(KQSP+LQFX+19)
      IQ(KQSP+LQFX+38) = IQ(KQSP+LQFX+20)
      GO TO 999

   45 IQUEST(11) = ISW
      CALL ZFATAM ('FZOPHR - error from IORITE.')
#endif
*
* fzophr44.inc
*
#if defined(CERNLIB_QMND3)
      IF (IACMOX.EQ.0)             GO TO 48
      NIL = 0
      ISW = MAGTP (1,LQ(LWR),LUNX,NWR,NIL)
      IF (ISW.EQ.0)                GO TO 999
      IQUEST(12) = MAGTP (20B,NIL,LUNX,NIL,NIL)
      IF (JBIT(IQUEST(12),10).EQ.0)          GO TO 45
      IF (IQ(KQSP+LQFX+37)+IQ(KQSP+LQFX+38).EQ.0)  GO TO 45
      IQ(KQSP+LQFX+37) = IQ(KQSP+LQFX+19)
      IQ(KQSP+LQFX+38) = IQ(KQSP+LQFX+20)
      GO TO 999

   45 IQUEST(11) = ISW
      CALL ZFATAM ('FZOPHR - error from MAGTP.')
#endif

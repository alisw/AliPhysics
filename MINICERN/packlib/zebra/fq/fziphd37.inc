*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*
* fziphd37.inc
*
#if defined(CERNLIB_QMIBMD)
      CALL ADREDKD (LUNI,LQ(LIN),NWMREC,JRECNO,ISW)
      IF (ISW.NE.0)                GO TO 840
#elif defined(CERNLIB_QMIRTD)
      CALL ADREDKD (LUNI,LQ(LIN),NWMREC,JRECNO,ISW)
      IF (ISW.NE.0)                GO TO 840
#elif 1
      READ (LUNI,REC=JRECNO,IOSTAT=ISW) (LQ(LIN+J-1),J=1,NWMREC)
      IF (ISW.NE.0)                GO TO 840
#endif

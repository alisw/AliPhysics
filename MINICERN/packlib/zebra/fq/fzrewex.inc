*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzrewex.inc
*
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
      IF (IACMOF.EQ.0)             GO TO 28
      CALL IORWND (LUNF,ISW)
      GO TO 29
#endif
*
* fzrewex.inc
*
#if defined(CERNLIB_QMND3)
      IF (IACMOF.EQ.0)             GO TO 28
      NIL = 0
      CALL MAGTP (13B,NIL,LUNF,NIL,NIL)
      GO TO 29
#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzclosex.inc
*
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
      IF (IACMOF.EQ.0)             GO TO 48
      CALL IOCLOS (LUNF,ISW)
      GO TO 49
#endif
*
* fzclosex.inc
*
#if defined(CERNLIB_QMND3)
      IF (IACMOF.NE.0)             GO TO 99
#endif

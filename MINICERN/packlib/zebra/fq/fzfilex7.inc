*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*
* fzfilex7.inc
*
#if defined(CERNLIB_QMCDC)
      IF (NEOFOP.LT.0)  NEOFOP = IOPTT + 2
*    --------------  sequences for output  ---------------------
#endif
*
* fzfilex7.inc
*
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
      IACMO = IOPTY
*    --------------  sequences for output  ---------------------
#endif
*
* fzfilex7.inc
*
#if defined(CERNLIB_QMND3)
      IF (IOPTT.NE.0)  THEN
          IF (IOPTF.EQ.0)  IOPTY = 1
          NEOF = 2
        ENDIF
      IACMO = IOPTY
*    --------------  sequences for output  ---------------------
#endif

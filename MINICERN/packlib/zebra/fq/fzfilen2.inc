*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*
* fzfilen2.inc
*
#if defined(CERNLIB_QMCDC)
      IF (NEOFOP.LT.0)  NEOFOP = IOPTT + 2
#endif
*
* fzfilen2.inc
*
#if defined(CERNLIB_QMND3)
      IF (IOPTT.EQ.1)  NEOF=2
#endif

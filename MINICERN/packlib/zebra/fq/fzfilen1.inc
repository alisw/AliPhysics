*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*
* fzfilen1.inc
*
#if defined(CERNLIB_QMLNX)
      IF (IOPTT.NE.0)  LREC=MIN(LREC,2499)
C-      maximum logical record size for tape files on the VAX (ANSI)
*    --------------  sequences for output  ---------------------
#endif
*    --------------  conversions for output  -------------------
*
* fzfilen1.inc
*
#if defined(CERNLIB_QMVAX)
      IF (IOPTT.NE.0)  LREC=MIN(LREC,2499)
C-      maximum logical record size for tape files on the VAX (ANSI)
*    --------------  sequences for output  ---------------------
#endif
*    --------------  conversions for output  -------------------
*
* fzfilen1.inc
*
#if defined(CERNLIB_QMVMI)
      IF (IOPTT.NE.0)  LREC=MIN(LREC,2499)
C-      maximum logical record size for tape files on the VAX (ANSI)
*    --------------  sequences for output  ---------------------
#endif
*    --------------  conversions for output  -------------------

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:54  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_QTRACE99_INC
#define CERNLIB_ZEBRA_QTRACE99_INC
*
*
* qtrace99.inc
*
#if !defined(CERNLIB_QTRHOLL)||defined(CERNLIB_A6M)
  999 NQTRAC = NQTRAC - 1
#else
  999 NQTRAC = NQTRAC - 2
#endif

#endif

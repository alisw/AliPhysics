*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*             Extra EOF test
*
* fzdaeofx.inc
*
#if defined(CERNLIB_QMIRT)
      IF (ISW.EQ.122)              GO TO 841
#endif
*             Extra EOF test
*
* fzdaeofx.inc
*
#if defined(CERNLIB_QMIRTD)
      IF (ISW.EQ.122)              GO TO 841
*    needs packing  machine words <-> file words
#ifndef CERNLIB_FQNEEDPK
#endif
*    needs convert  machine words <-> file words
#ifndef CERNLIB_FQNEEDCV
#endif
*    --------------  sequences for output  ---------------------
#endif

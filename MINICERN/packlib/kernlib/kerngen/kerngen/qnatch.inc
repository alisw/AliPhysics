*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:14  mclareni
* Kernlib
*
*
*  Siemens compiler
*
* qnatch.inc
*
#if defined(CERNLIB_QMIBM)
#if defined(CERNLIB_QF_SIEM)
      COMMON /QNATCH/ NATCH(256)
      INTEGER*2       NATCH
#else
      INTEGER*2       NATCH(256)
#endif
#include "qnatchdata.inc"
#endif

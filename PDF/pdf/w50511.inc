*
* $Id$
*
* $Log$
* Revision 1.4  2000/09/18 10:02:34  hristov
* Makefile added to PDF8
*
* Revision 1.1.1.2  1996/10/30 08:27:12  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:49  plothow
* Version 7.01
*
*
*
* w50511.inc
*
      COMMON/W50511/ NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION
#endif
#if defined(CERNLIB_SINGLE)
      REAL
#endif
     +       TMAS

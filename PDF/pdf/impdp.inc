*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:27:10  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:48  plothow
* Version 7.01
*
*
*
* impdp.inc
*
#if defined(CERNLIB_DOUBLE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#endif
#if defined(CERNLIB_SINGLE)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#endif

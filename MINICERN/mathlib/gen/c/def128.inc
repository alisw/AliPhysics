*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:12  mclareni
* Mathlib gen
*
*
*
* def128.inc
*
#if (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
      REAL(4)
#endif
#if (!defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
      REAL(2)
#endif
#if (!defined(CERNLIB_F90))&&(defined(CERNLIB_NOQUAD)||!defined(CERNLIB_DOUBLE))
      DOUBLE PRECISION
#endif
      REAL*16

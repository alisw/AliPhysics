*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp128.inc
*
#if (!defined(CERNLIB_F90))&&(defined(CERNLIB_NOQUAD)||!defined(CERNLIB_DOUBLE))
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#elif (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      IMPLICIT REAL(2) (A-H,O-Z)

#elif (defined(CERNLIB_F90))&&(defined(CERNLIB_DOUBLE))
      IMPLICIT REAL(4) (A-H,O-Z)

#elif 1
      IMPLICIT REAL*16 (A-H,O-Z)

#endif

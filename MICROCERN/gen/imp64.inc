*
* $Id$
*
* $Log$
* Revision 1.1.2.1  2002/07/11 17:14:49  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:36  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
#if !defined(CERNLIB_DOUBLE)
      IMPLICIT REAL (A-H,O-Z)

#elif (defined(CERNLIB_F90))&&(defined(CERNLIB_DOUBLE))
      IMPLICIT REAL(2) (A-H,O-Z)

#elif 1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#endif

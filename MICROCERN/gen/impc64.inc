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
*
* impc64.inc
*
#if (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
      IMPLICIT COMPLEX(2) (A-H,I,O-Z)

#elif !defined(CERNLIB_DOUBLE)
      IMPLICIT COMPLEX    (A-H,I,O-Z)

#elif 1
      IMPLICIT COMPLEX*16 (A-H,I,O-Z)

#endif

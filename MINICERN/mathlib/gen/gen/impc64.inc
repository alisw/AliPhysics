*
* $Id$
*
* $Log$
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

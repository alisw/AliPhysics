*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* defc64.inc
*
#if (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
      COMPLEX(2)

#elif !defined(CERNLIB_DOUBLE)
      COMPLEX

#elif 1
      COMPLEX*16

#endif

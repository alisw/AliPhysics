*
* $Id$
*
* $Log$
* Revision 1.1.2.1  2002/07/11 17:14:48  alibrary
* Adding MICROCERN
*
* Revision 1.1.1.1  1999/05/18 15:55:35  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*
#if !defined(CERNLIB_DOUBLE)
      REAL
#endif
#if (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
      REAL(2)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION
#endif

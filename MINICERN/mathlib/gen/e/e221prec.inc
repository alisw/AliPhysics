*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
*
* e221prec.inc
*
#if defined(CERNLIB_IBM)
      DATA BIG /1.E+75/

#elif (defined(CERNLIB_SINGLE))&&(defined(CERNLIB_CONVEX))
      DATA BIG /1.E+307/

#elif 1
      DATA BIG /1.E+37/

#endif

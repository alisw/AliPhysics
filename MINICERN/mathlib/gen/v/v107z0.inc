*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/02 22:53:48  thakulin
* Add support for Apogee and EPC Fortran compilers.
*
* Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
* Mathlib gen
*
*
*
* v107z0.inc
*

#if defined(CERNLIB_ALLIANT)||defined(CERNLIB_DECS)||defined(CERNLIB_MACMPW)||defined(CERNLIB_SGI)
      DATA IZ0/'33000000'X/

#elif defined(CERNLIB_APOLLO)
      DATA IZ0/16#33000000/

#elif defined(CERNLIB_IBMRT)||defined(CERNLIB_IBMAIX)
      DATA IZ0/Z33000000/

#elif (defined(CERNLIB_CONVEX)||defined(CERNLIB_HPUX)||defined(CERNLIB_SUN))&&!(defined(CERNLIB_F90)||defined(CERNLIB_QFEPC))
      DATA IZ0/X'33000000'/

#elif 1
      DATA IZ0/Z'33000000'/

#endif

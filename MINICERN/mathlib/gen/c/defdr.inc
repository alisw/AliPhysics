*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:12  mclareni
* Mathlib gen
*
*
*
* defdr.inc
*
#if (defined(CERNLIB_QF2C))&&(!defined(CERNLIB_LINUX))
      DOUBLE PRECISION DREAL,DIMAG,DCONJG,DCMPLX
      EXTERNAL         DREAL,DIMAG,DCONJG,DCMPLX
#endif
#if (defined(CERNLIB_QF2C))&&(defined(CERNLIB_LINUX))
      DOUBLE PRECISION DREAL
      EXTERNAL         DREAL
#endif

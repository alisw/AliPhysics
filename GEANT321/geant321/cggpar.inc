*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGGPAR_INC
#define CERNLIB_GEANT321_CGGPAR_INC
*
*
* cggpar.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (LCGHEA=3, KCGSIZ=1)
      PARAMETER (          KCGATT=2)
      PARAMETER (          KCGNF =3)
      PARAMETER (LCGFAC=6, KCGAF =1)
      PARAMETER (          KCGAA =2)
      PARAMETER (          KCGBB =3)
      PARAMETER (          KCGCC =4)
      PARAMETER (          KCGDD =5)
      PARAMETER (          KCGNE =6)
      PARAMETER (LCGEDG=7, KCGAE =1)
      PARAMETER (          KCGX1 =2)
      PARAMETER (          KCGY1 =3)
      PARAMETER (          KCGZ1 =4)
      PARAMETER (          KCGX2 =5)
      PARAMETER (          KCGY2 =6)
      PARAMETER (          KCGZ2 =7)
#endif
*
* cggpar.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGGPAR.
#include "geant321/incend.inc"
*
#endif

#endif

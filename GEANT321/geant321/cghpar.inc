*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGHPAR_INC
#define CERNLIB_GEANT321_CGHPAR_INC
*
*
* cghpar.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (I4SIZE = 1)
      PARAMETER (LHHEAD=18,KHSIZE=1)
      PARAMETER (          KHNT  =2)
      PARAMETER (          KHIFWI=3)
      PARAMETER (          KHNXYZ=8)
      PARAMETER (          KHNFAC=9)
      PARAMETER (          KHNEDG=10)
      PARAMETER (          KHJXYZ=11)
      PARAMETER (          KHJFA =12)
      PARAMETER (          KHJPFA=13)
      PARAMETER (          KHJDFA=14)
      PARAMETER (          KHJTRE=15)
      PARAMETER (          KHJSTA=16)
      PARAMETER (          KHJALE=17)
      PARAMETER (          KHJARI=18)
#endif
*
* cghpar.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGHPAR.
#include "geant321/incend.inc"
*
#endif

#endif

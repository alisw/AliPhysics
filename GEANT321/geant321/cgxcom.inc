*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGXCOM_INC
#define CERNLIB_GEANT321_CGXCOM_INC
*
*
* cgxcom.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON /CGXCOM/ IXLUN,NXMAX,NXCUR,NXNODE,NXEDGE,NXFACE
#endif
*
* cgxcom.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGXCOM.
#include "geant321/incend.inc"
*
#endif

#endif

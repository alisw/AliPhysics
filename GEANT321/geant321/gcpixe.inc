*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:36  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCPIXE_INC
#define CERNLIB_GEANT321_GCPIXE_INC
*
*
* gcpixe.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCPIXE/LIMPRE,IFLAPE,ICOLOR,IXXX,IYYY,
     +ISSEEN,ISCOLO,ISLSTY,ISLWID,ISFILL,
     +IMAP,JON,NMAP,UUU,VVV,ZUV,ZNMAP1
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtpixe.inc"
#endif
 
#endif

#endif

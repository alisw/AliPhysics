*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:34  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCPARA_INC
#define CERNLIB_GEANT321_GCPARA_INC
*
*
* gcpara.inc
*
#if !defined(CERNLIB_INCLUDE)
#include "geant321/gtpara.inc"
      PARAMETER (LSTACK = 5000)
C     BITPOT is for Phi.Or.Tet
C
C ---------------------------------------------------------
      COMMON    /GCPARA/
     +                   EPSIX0 (LSTACK)       ,
     +                   IDRPHI (LSTACK     )  , IDRTET (LSTACK     ),
     +                   IDROUT (LSTACK     )  , JPLOST (LSTACK     ),
     +                   IPHTMP (LSTACK     )  ,
     +                   BITPHI (LSTACK     )  , BITTET (LSTACK     ),
     +                   BITPOT (LSTACK     )  , JJLOST, JJFILL,
     +                                           JENTRY, JEMPTY,
     +                                           EPSMAX,
     +                   JJTEMP, JJWORK        , JJSTK1,
     +                   J1TEMP,                 J1STK1,
     +                   IFOUNP, IFOUNT        , IFNPOT,
     +                                           SYMPHI,
     +                   SYMTEU, SYMTED
C
#endif
*
* gcpara.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcpara
#include "geant321/incend.inc"
#endif

#endif

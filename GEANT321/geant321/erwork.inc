*
* $Id$
*
* $Log$
* Revision 1.2  1996/07/25 15:01:16  ravndal
* CMZ to CVS conversion artefact removed
*
* Revision 1.1  1996/03/06 15:32:34  mclareni
* Add geane321 includes
*
*
#ifndef CERNLIB_ERPRCDES_ERWORK_INC
#define CERNLIB_ERPRCDES_ERWORK_INC
*
*      Working sequences for the interface
*
* erwork.inc
*
      DOUBLE PRECISION EI, EF, ASDSC
      COMMON /ERWORK/ EI(15), EF(15), ASDSC(5,5),
     +                   XI(3), PPI(3), HI(9),
     +                   XF(3), PF(3),  HF(9),
     +                   CHTR, DEDX2, BACKTR, CUTEK, TLGCM2, TLRAD
*
* erwork.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
     +    ERWORK
#include "geant321/incend.inc"
#endif

#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_AUXPAR_INC
#define CERNLIB_GEANT321_AUXPAR_INC
*
*
* auxpar.inc
*
C***********************************************************************
C     INCLUDE FILE: AUXPAR COPY                   CREATED 11/12/86 BY PA
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX6=1133)
      CHARACTER*8 ANA
      COMMON /FKAUXP/PXA(IDMAX6),PYA(IDMAX6),PZA(IDMAX6),HEPA(IDMAX6),
     *AMA(IDMAX6),ICHA(IDMAX6),IBARA(IDMAX6),            NREA(IDMAX6)
      COMMON /FKAUXC/ANA(IDMAX6)
 

#endif

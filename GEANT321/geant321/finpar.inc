*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:39  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_FINPAR_INC
#define CERNLIB_GEANT321_FINPAR_INC
*
*
* finpar.inc
*
C***********************************************************************
C     INCLUDE FILE: FINPAR COPY                   CREATED 11/12/86 BY PA
C     NOTE: DO NOT FORGET FILE FINPAR2 COPY
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX1=1133)
      CHARACTER*8 ANF
      COMMON /FKFINP/PXF(IDMAX1),PYF(IDMAX1),PZF(IDMAX1),HEP(IDMAX1),
     *AMF(IDMAX1),ICHF(IDMAX1),IBARF(IDMAX1),NREF(IDMAX1)
      COMMON /FKFINC/ ANF(IDMAX1)
 

#endif

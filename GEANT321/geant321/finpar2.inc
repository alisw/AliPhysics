*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:39  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_FINPAR2_INC
#define CERNLIB_GEANT321_FINPAR2_INC
*
*
* finpar2.inc
*
C***********************************************************************
C     INCLUDE FILE: FINPAR2 COPY                   CREATED 11/12/86 BY A
C     NOTE: DO NOT FORGET FILE FINPAR COPY
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX3=1133)
      CHARACTER*8 ANF
      COMMON /FKFINP/PXF(IDMAX3),PYF(IDMAX3),PZF(IDMAX3),HEF(IDMAX3),
     *AMF(IDMAX3),ICHF(IDMAX3),IBARF(IDMAX3),NREF(IDMAX3)
      COMMON /FKFINC/ ANF(IDMAX3)
 

#endif

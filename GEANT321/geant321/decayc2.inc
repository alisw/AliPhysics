*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:38  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_DECAYC2_INC
#define CERNLIB_GEANT321_DECAYC2_INC
*
*
* decayc2.inc
*
C***********************************************************************
C     INCLUDE FILE: DECAYC2 COPY                   CREATED 12/12/86 BY A
C     NOTE: SEE ALSO DECAYC COPY
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX9=540)
      CHARACTER*8 ZKNAMC
      COMMON /FKDECC/ NZKC(IDMAX9,3),WTC(IDMAX9)
      COMMON/FKDECH/ ZKNAMC(IDMAX9)
 

#endif

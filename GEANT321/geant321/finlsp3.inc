*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:38  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_FINLSP3_INC
#define CERNLIB_GEANT321_FINLSP3_INC
*
*
* finlsp3.inc
*
C***********************************************************************
C     include file:  finlsp3 copy           Created 4/9/89 by A. Ferrari
C     note: see also file finlsp copy and finlsp2
C     changes: none
C     included in the following subroutines or functions:
C                                                         Eventq
C
C     description of the common block(s) and variable(s)
C
C***********************************************************************
      PARAMETER (IDMAX7=100)
      COMMON/FKFINL/CXH(IDMAX7),CYH(IDMAX7),CZH(IDMAX7),
     *ELH(IDMAX7),PLH(IDMAX7),IH,ITH(IDMAX7)
 

#endif

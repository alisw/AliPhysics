*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      REAL FUNCTION FCTCOS(T)
C
C *** NVE 01-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (27-OCT-1983)
C
#include "geant321/s_coscom.inc"
#include "geant321/limits.inc"
C
      DOUBLE PRECISION TEST1,TEST2
C
      TEST1=-BB*T*1.0D0
      IF (TEST1 .GT. EXPXU) TEST1=EXPXU
      IF (TEST1 .LT. EXPXL) TEST1=EXPXL
      TEST2=-DD*T*1.0D0
      IF (TEST2 .GT. EXPXU) TEST2=EXPXU
      IF (TEST2 .LT. EXPXL) TEST2=EXPXL
C
      FCTCOS=AA*EXP(TEST1)+CC*EXP(TEST2)-RR
C
      RETURN
      END

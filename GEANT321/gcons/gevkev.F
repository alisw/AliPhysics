*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GEVKEV(EGEV,ENERU,KUNIT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Convert energy from GeV to more appropriate units        *
C.    *                                                                *
C.    *    ==>Called by : GPLMAT,GPRMAT,GPHYSI,GPCXYZ                  *
C.    *       Author     M.Maire    *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
      CHARACTER*4 KUNIT
*
      IF(EGEV.LT.0.0000009999) THEN
         ENERU = EGEV*1.E+9
         KUNIT = '  eV'
      ELSEIF(EGEV.LT.0.0009999) THEN
         ENERU = EGEV*1.E+6
         KUNIT = ' keV'
      ELSEIF(EGEV.LT.0.9999) THEN
         ENERU = EGEV*1.E+3
         KUNIT = ' MeV'
      ELSEIF(EGEV.LT.999.9) THEN
         ENERU = EGEV
         KUNIT = ' GeV'
      ELSE
         ENERU = EGEV*1.E-3
         KUNIT = ' TeV'
      ENDIF
*
      END

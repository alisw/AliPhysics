*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GIDROP
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      DROP all ZEBRA structures created at initialisation time  *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF(JMATE .NE.0)CALL MZDROP(IXCONS,JMATE,' ')
      IF(JTMED .NE.0)CALL MZDROP(IXCONS,JTMED,' ')
      IF(JVOLUM.NE.0)CALL MZDROP(IXCONS,JVOLUM,'L')
      IF(JPART .NE.0)CALL MZDROP(IXCONS,JPART,' ')
      IF(JROTM .NE.0)CALL MZDROP(IXCONS,JROTM,' ')
      IF(JSET  .NE.0)CALL MZDROP(IXCONS,JSET,' ')
      JVOLUM=0
      JMATE=0
      JTMED=0
      JPART=0
      JROTM=0
      JSET=0
*
      CALL VZERO (NLEVEL,306)
*
      END
 

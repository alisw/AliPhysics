*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:40  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GDEBUG
*
************************************************************************
*             Standard GEANT debug                                     *
*                                                                      *
*          ==>Called by : <USER>,GUSTEP                                *
*             Author    R.Brun  *********                              *
************************************************************************
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
*
*             Debug/plot event
      IF(IDEBUG.NE.0) THEN
        IF((ISWIT(2).EQ.1).OR.(ISWIT(3).EQ.1)
     +                    .OR.(ISWIT(2).EQ.4)) CALL GSXYZ
        IF (ISWIT(2).EQ.2) CALL GPCXYZ
        IF (ISWIT(1).EQ.2) CALL GPGKIN
        IF (ISWIT(2).EQ.3) THEN
           IF(ISWIT(4).EQ.3.AND.CHARGE.EQ.0.)RETURN
           CALL GDCXYZ
           CALL IGSA(0)
        ELSEIF(ISWIT(2).EQ.4) THEN
           IF(ISTOP.NE.0) CALL GDTRAK(' ')
        ENDIF
      ENDIF
      END
 

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE FRADNC.FOR
*COPY FRADNC
*
*=== fradnc ===========================================================*
*
      FUNCTION FRADNC (RHO)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
*
#include "geant321/nucstf.inc"
*
      IF ( RHO .GE. RHOCOR ) THEN
         FRADNC = RACORE (RHO)
      ELSE IF ( RHO .GE. RHOSKN ) THEN
         FRADNC = RASKIN (RHO)
      ELSE IF ( RHO .GT. 0.D+00 ) THEN
         FRADNC = RAHALO (RHO)
      ELSE
         FRADNC = RADTOT
      END IF
      RETURN
      END

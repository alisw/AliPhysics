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
*$ CREATE FPFRNC.FOR
*COPY FPFRNC
*
*=== fpfrnc ===========================================================*
*
      FUNCTION FPFRNC (RHO,I)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
*
#include "geant321/nucstf.inc"
*
      IF ( RHO .GE. RHOCOR ) THEN
         FPFRNC = PFCORE (RHO,I)
      ELSE IF ( RHO .GE. RHOSKN ) THEN
         FPFRNC = PFSKIN (RHO,I)
      ELSE IF ( RHO .GT. 0.D+00 ) THEN
         FPFRNC = PFHALO (RHO,I)
      ELSE
         FPFRNC = 0.D+00
      END IF
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE DRELAB.FOR
*COPY DRELAB
*
*=== drelab ===========================================================*
*
      SUBROUTINE DRELAB(X,Y,Z,COTE,SITE,COPS,SIPS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      X1=COPS*X-SIPS*COTE*Y+SIPS*SITE*Z
      X2=SIPS*X+COPS*COTE*Y-COPS*SITE*Z
      X3=SITE*Y+COTE*Z
      X=X1
      Y=X2
      Z=X3
      RETURN
      END

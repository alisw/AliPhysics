*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE SFECFE.FOR
*COPY SFECFE
*
*=== sfecfe ===========================================================*
*
      SUBROUTINE SFECFE(SFE,CFE)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      REAL RNDM(2)
      ENTRY         COSI(SFE,CFE)
C********************************************************************
C     VERSION JUNE 81 BY             PERTTI AARNIO
C     LAST CHANGE 11. DECEMBER 85 BY PERTTI AARNIO
C                                    HELSINKI UNIVERSITY OF
C                                    TECHNOLOGY, FINLAND
C
C
C     SUBROUTINE OF FLUKA TO GIVE SINE AND COSINE OF AN
C     RANDOM ANGLE UNIFORMLY DISTRIBUTED BETWEEN 0 AND 2*PI
C********************************************************************
C
 10   CALL GRNDM(RNDM,2)
      X=2.D0*RNDM(1)-1.D0
      Y=RNDM(2)
      X2=X*X
      Y2=Y*Y
      IF (X2+Y2.GT.1.D0) GO TO 10
      CFE=(X2-Y2)/(X2+Y2)
      SFE=2.D0*X*Y/(X2+Y2)
      RETURN
      END

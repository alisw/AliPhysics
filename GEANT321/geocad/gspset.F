*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GSPSET(PAR)
C-
C-   Created  12-JUN-1991   Nils Joar Hoimyr
C-   Modified 21.02.1992  Jouko Vuoskoski
C-
C-  Takes a GEANT SPHE shape, and wites it out as a primitive solid
C-  sphere in SET.
C----------------------------------------------------------------------
 
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(50)
C
      REAL RMIN, RMAX
C
C Segments not yet ready:
C    REAL ALPH1, ALPH2, THET1, THET2
C
C----------------------------------------------------------------------
C
      RMIN= PAR(1)
      RMAX= PAR(2)
C
C        IF segment of a sphere  THEN (Currently not implemented)
C      *WRITE SET @50,N..#34,RMAX,ALPH1,ALPH2,THET1,THET2,RMIN
C        ELSE
C      *WRITE SET @50,N..#34,RMAX,,,,,RMIN
C----------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1,RMAX,RMIN
      CALL GJWRIT
C
10000   FORMAT('@50,',I10,',:5,2#34,',G14.7,',,,,,',G14.7)
C
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GBOSET(PAR)
C-   Created  12-JUN-1991   Nils Joar Hoimyr, Jouko Vuoskoski
C-   Modified 21.02.1992  Jouko Vuoskoski
C-
C- Takes the GEANT shape BOX, converts it to a SET box, with
C- the right position and writes the shape to the current SET file.
C -
C----------------------------------------------------------------------
C
#include "geant321/gcsetf.inc"
 
C
      REAL      PAR
      DIMENSION PAR(50)
C
      REAL X,Y,Z,DX,DY,DZ
C
C
C----------------------------------------------------------------------
C
C             BOX
C
      DX=PAR(1)
      DY=PAR(2)
      DZ=PAR(3)
      X=2*DX
      Y=2*DY
      Z=2*DZ
C *WRITE SET @50,N1,:5,2#60,X,Y,Z
C *WRITE SET @302,N2,#317,-DX,-DY,-DZ
C *WRITE SET @100,N3,:5,2,:9,'Name'#101,!N1,!N2
C
C-------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1, X, Y, Z
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10100)N1,-DX,-DY,-DZ
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10200)N1,N1-2,N1-1
      CALL GJWRIT
C
10000   FORMAT('@50,',I10,',:5,2#60,',G14.7,',',G14.7,',',G14.7)
10100   FORMAT('@302,',I10,'#317,',G14.7,',',G14.7,',',G14.7)
10200   FORMAT('@100,',I10,',:5,2#101,!',I10,',!',I10)
C
C
C
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GELSET(PAR)
C-
C-   Created  26-JUL-1991   Nils Joar Hoimyr
C-   Modified 21.02.1992  Jouko Vuoskoski
C-
C-   Writes out the elliptical section of a tube in SET, and generates
C-   an extruded solid from the elliptical section (face).
C-
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(100)
      REAL DX1,DY,Z,DZ,A,B
C
C------------------------------------------------------
C  Converts GEANT shape ELTU into SET format
C
      A=PAR(1)
      B=PAR(2)
      DZ=PAR(3)
C       SET conversion
      Z= 2*DZ
      DX1=0.0
      DY=0.0
C *WRITE SET @103,N1..#11,A,B,0.0,360.0,0,0,0,0
C *WRITE SET @100,N2,:5,2#143,!N!,0.0,0.0,1.0
C *WRITE SET @302,N3,:5,2#317,0,0,-DZ
C *WRITE SET @100,N4,:5,2,:9,'MATNAM'#101,!N1,!N2
C
C------------------------------------------------
C
      WRITE(BLKSTR,10000)N1,A,B
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10100)N1,N1-1,Z
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10200)N1,-DX1,-DY,-DZ
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10300)N1,N1-2,N1-1
      CALL GJWRIT
10000   FORMAT('@103,',I10,',:5,2#11,:14,1,',G14.7,',',G14.7,',
     +     0.0,360.0,0.0,0.0,0.0,0.0')
10100   FORMAT('@100,',I10,',:5,2#143,!',I10,',0.,0.,',G14.7)
10200   FORMAT('@302,',I10,'#317,',G14.7,',',G14.7,',',G14.7)
10300   FORMAT('@100,',I10,',:5,2#101,!',I10,',!',I10)
C
   10 RETURN
      END

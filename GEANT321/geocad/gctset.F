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
      SUBROUTINE GCTSET(PAR)
C-
C-   Created  26-JUL-1991   Nils Joar Hoimyr
C-   Modified by J. Vuoskoski 21.02.1992
C-
C-   Describes a contour face of a tube segment from the GEANT
C-   CTUB shape parameters.  This face is rotated around the
C-   Z-axis to generated a solid tube segment in SET.  To cut the
C-   tube, the cutting planes are calculated from the shape parameters,
C-   and then 2 half-spaces are generated and subtracted from the tube
C-   in a boolean CUT operation.
C
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(100)
C
      REAL  Z,DX,DY,DZ,RMIN,RMAX,PHIMIN,PHIMAX
      REAL  LXL,LYL,LZL,LXH,LYH,LZH
C
C----------------------------------------------------------
C
      RMIN=PAR(1)
      RMAX=PAR(2)
      DX=0.0
      DY=0.0
      DZ=PAR(3)
      PHIMIN=PAR(4)
      PHIMAX=PAR(5)
      LXL= PAR(6)
      LYL= PAR(7)
      LZL= PAR(8)
      LXH= PAR(9)
      LYH= PAR(10)
      LZH= PAR(11)
      IF (RMAX .LE. 0.0) THEN
         WRITE (*,*) 'IMPOSSIBLE RADIUS VALUE'
         N1=N1-1
         GOTO 10
      ENDIF
C
C       SET CONVERSION
      Z= 2*DZ
C
C   Starts with a normal TUBS element
C *WRITE SET @50,N1,:5,2#32,RMAX,Z,PHIMIN,PHIMAX,RMIN
C *WRITE SET @302,N2,:5,2#317,0,0,-DZ
C *WRITE SET @100,N3,:5,2,:9,'MATNAM'#101,!N1,!N2
C
C------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1,RMAX,Z,PHIMIN,PHIMAX,RMIN
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10100)N1,-DX,-DY,-DZ
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10200)N1,N1-2,N1-1
      CALL GJWRIT
C
C     Tube cutting:
C
C     Create cutting planes, and semi-spaces which are used to
C     cut the tube section
C     SET blocks @30..#30 plane normals given by LXL, etc
C
      N1=N1+1
      WRITE(BLKSTR,10300)N1,LXL,LYL,LZL,-DZ
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10400)N1,N1-1
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10300)N1,LXH,LYH,LZH,DZ
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10400)N1,N1-1
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10500)N1,N1-5,N1-1,N1-3
      CALL GJWRIT
C
10000   FORMAT('@50,',I10,',:5,2#32,',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7,',',G14.7)
10100   FORMAT('@302,',I10,'#317,',G14.7,',',G14.7,',',G14.7)
10200   FORMAT('@100,',I10,',:5,2#101,!',I10,',!',I10)
10300   FORMAT('@30,',I10,',:5,2#30,',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7)
10400   FORMAT('@100,',I10,',:5,2#139,!',I10,',1,1')
10500   FORMAT('@100,',I10,',:5,2#100,3,!',I10,',!',I10,',!',I10)
C
   10 RETURN
      END

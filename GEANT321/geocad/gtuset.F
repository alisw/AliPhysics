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
      SUBROUTINE GTUSET(PAR,ISHAPE)
C-
C-   Created  12-JUN-1991   Nils Joar Hoimyr
C-   Modified 21.02.1992  Jouko Vuoskoski
C-
C   Converts GEANT shapes TUBE&TUBS to SET
C-
C-   Describes a contour face of a TUBE(Segment) from the GEANT
C-   TUBE/TUBS shape parameters.  This face is rotated around the
C-   Z-axis to generated the solid shape in SET.
C
C---------------------------------------------------------------------
C
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(50)
      REAL DX1,DZ,RMIN,RMAX,PHIMIN,PHIMAX,Z2
C
C----------------------------------------------------------------------
C
      RMIN=PAR(1)
      RMAX=PAR(2)
      Z2=PAR(3)
      IF (ISHAPE.EQ.5) THEN
         PHIMIN=0.0
         PHIMAX=360.0
      ELSEIF (ISHAPE.EQ.6) THEN
         PHIMIN=PAR(4)
         PHIMAX=PAR(5)
      ENDIF
C       SET CONVERSION
      Z= 2*Z2
      DX1=0.0
      DY=0.0
      DZ=Z2
      IF (RMAX .LE. 0.0) THEN
         WRITE (*,*) 'IMPOSSIBLE RADIUS VALUE'
         N1=N1-1
         GOTO 10
      ENDIF
C *WRITE SET @50,N1,:5,2#32,RMAX,Z,PHIMIN,PHIMAX,RMIN
C *WRITE SET @302,N2,:5,2#317,0,0,-Z2
C *WRITE SET @100,N3,:5,2,:9,'MATNAM'#101,!N1,!N2
C-------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1,RMAX,Z,PHIMIN,PHIMAX,RMIN
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10100)N1,-DX1,-DY,-DZ
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10200)N1,N1-2,N1-1
      CALL GJWRIT
C
10000   FORMAT('@50,',I10,',:5,2#32,',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7,',',G14.7)
10100   FORMAT('@302,',I10,'#317,',G14.7,',',G14.7,',',G14.7)
10200   FORMAT('@100,',I10,',:5,2#101,!',I10,',!',I10)
C
   10 RETURN
      END

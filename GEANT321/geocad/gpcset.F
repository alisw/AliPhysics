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
      SUBROUTINE GPCSET(PAR)
C-
C-   Created  26-JUN-1991   Nils Joar Hoimyr
C-   Modified 21.02.1992  Jouko Vuoskoski
C-   Modified  02.04.92  NJH
C-
C---------------------------------------------------------
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(100)
      REAL PHIMIN,PHIMAX,PX,PY,PZ,RMIN,RMAX
C-
C-   Describes a contour face of a PCON shape calculated from the GEANT
C-   PCON shape parameters.  This face is rotated around the
C-   Z-axis to generated the solid shape in SET.
C---------------------------------------------------------
C
      PHIMIN=(PAR(1)-90.0)
      PHIMAX=(PHIMIN+PAR(2))
      NZ=PAR(3)
C        Solid of revolution described by rotation of a face
C        Face defined in the yz-plane (x=0)
C       WRITE SET@103,F1#3,3,2,P1X,P1Y,P1Z,P2X.........
C---------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1
      CALL GJWRIT
C---------------------------------------------------------
      N2= 1
      PX= 0.0
C
      PY= PAR(5)
      PZ= PAR(4)
C---------------------------------------------------------
C
      WRITE(BLKSTR,10100)PX,PY,PZ
      CALL GJWRIT
C---------------------------------------------------------
C
      DO 10  K=1,NZ
         N2=N2+3
         PY= PAR(N2+2)
         PZ= PAR(N2)
         RMIN= PAR(N2+1)
         RMAX=PAR(N2+2)
         IF (RMIN .GE. RMAX) GOTO 10
C
C---------------------------------------------------------
C
         WRITE(BLKSTR,10100)PX,PY,PZ
         CALL GJWRIT
C---------------------------------------------------------
C
   10 CONTINUE
C
      DO 20  L=2,NZ
         PY= PAR(N2+1)
         PZ= PAR(N2)
C
C---------------------------------------------------------
C
         WRITE(BLKSTR,10100)PX,PY,PZ
         CALL GJWRIT
C---------------------------------------------------------
         N2=N2-3
   20 CONTINUE
C     WRITE SET @100,N..#144,!F1,PHIMIN,PHIMAX.
C
      N1=N1+1
      WRITE(BLKSTR,10200)N1,N1-1,PHIMIN,PHIMAX
      CALL GJWRIT
C
10000   FORMAT('@103,',I10,',:5,2#3,3,2')
10100   FORMAT(',',G14.7,',',G14.7,',',G14.7)
10200   FORMAT('@100,',I10,',:5,2#144,!',I10,',',G14.7,',',G14.7)
C
      END

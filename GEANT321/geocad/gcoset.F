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
      SUBROUTINE GCOSET(PAR,ISHAPE)
C-
C-   Created  12-JUN-1991   Nils Joar Hoimyr
C-   Modified 21.02.1992    Jouko Vuoskoski
C-
C-   Describes a contour face of a CONE(Segment) from the GEANT
C-   CONE/CONS shape parameters.  This face is rotated around the
C-   Z-axis to generated the solid shape in SET.
C----------------------------------------------------------C
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(100)
      REAL P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,P3Y,P3Z,P4X,P4Y,P4Z
C
C---------------------------------------------------------
      RMIN1=PAR(2)
      RMAX1=PAR(3)
      RMIN2=PAR(4)
      RMAX2=PAR(5)
      Z2=PAR(1)
      Z1=-Z2
      IF (ISHAPE.EQ.7) THEN
         PHIMIN=0.0
         PHIMAX=360.0
      ELSEIF (ISHAPE.EQ.8) THEN
         PHIMIN=PAR(6)
         PHIMAX=PAR(7)
      ENDIF
C       SET CONVERSION
      P1X= RMIN1
      P1Y= 0
      P1Z= Z1
      P2X= RMIN2
      P2Y= 0
      P2Z= Z2
      P3X= RMAX2
      P3Y= 0
      P3Z= Z2
      P4X= RMAX1
      P4Y= 0
      P4Z= Z1
C *WRITE SET@103,F1#3,3,2,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,P3Y,P3Z,P4X,P4Y,P4Z
C *WRITE SET @100,N..#144,!F1,PHIMIN,PHIMAX.
C----------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1,P1X,P1Y,P1Z,P2X,P2Y,P2Z, P3X,P3Y,P3Z,P4X,
     +P4Y,P4Z
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10100)N1,N1-1,PHIMIN,PHIMAX
      CALL GJWRIT
C
10000   FORMAT('@103,',I10,',:5,2#3,3,2,',G14.7,',',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7,',',G14.7)
10100   FORMAT('@100,',I10,',:5,2#144,!',I10,',',G14.7,',',G14.7)
C
C
      END

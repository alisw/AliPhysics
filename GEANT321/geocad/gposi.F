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
      SUBROUTINE GPOSI(PARMAT,POSIT,VOLNAM,COLO)
C-   Created  12-JUN-1991   Nils Joar Hoimyr, Jouko Vuoskoski
C
C    Modyfied 10.11.92 Nils Hoimyr. Added colour output to CAD.
C----------------------------------------------------------------------
C
#include "geant321/gcsetf.inc"
 
C
C
      DIMENSION PARMAT(9)
      DIMENSION POSIT(3)
      CHARACTER*10 VOLNAM
      INTEGER COLOUR,COL,COLO
C
C
C    Map GEANT colour to SET colour: (Use SET standard colour table)
C                                     Approx LightBlue with blue
C
      IF (ABS(COLO).GE.10) THEN
         COL= ABS(COLO/10)
      ELSE
         COL= ABS(COLO)
      ENDIF
      IF (COL.EQ.1) COLOUR=8
      IF (COL.EQ.2) COLOUR=7
      IF (COL.EQ.3) COLOUR=4
      IF (COL.EQ.4) COLOUR=3
      IF (COL.EQ.5) COLOUR=5
      IF (COL.EQ.6) COLOUR=2
      IF (COL.EQ.7) COLOUR=3
C----------------------------------------------------------------------
C   Positioning of the volume in SET:
C
C *WRITE SET @302,N1,#301,1,0,0,0,1,0,0,0,1,XPOSIT,YPOSIT,XPOSIT
C *WRITE SET @100,N2,:57,'Colour',:9,'Name'#101,!N1(of the shape),!N1
C
C-------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1,(PARMAT(JV),JV=1,9),(POSIT(JV),JV=1,3)
      CALL GJWRIT
      N1=N1+1
      WRITE(BLKSTR,10100)N1,COLOUR,VOLNAM,N1-2,N1-1
      CALL GJWRIT
C
10000   FORMAT('@302,',I10,'#301,',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7,','
     +  ,G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7)
C10010   FORMAT('@302,',I10,'#317,',G14.7,',',G14.7,',',G14.7)
10100   FORMAT('@100,',I10,',:57,',I10,',:9,'''
     +  ,A,'''#101,!',I10,',!',I10)
C10020   FORMAT('@100,',I10,',:5,2#101,!',I10,',!',I10)
C
      END

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
      SUBROUTINE GPGSET(PAR)
C-
C-   Created  26-JUN-1991   Nils Joar Hoimyr
C-   Modified 21.02.1992  Jouko Vuoskoski
C-
C---------------------------------------------------------
C- Calculates the sectional face of 1 cell of the PGON shape.   This face is
C- copied and rotated around the Z-axis to make "the other side" of the first
C- cell.  A ruled solid is generated between the 2 faces, to created the cell.
C- The rest of the cells are created by copying and rotation of the first cell.
C- The final result is obtained with a Boolean fusion of the cells.  (The whole
C-  creation history is written to the SET file.)
C----------------------------------------------------------
C
#include "geant321/gcsetf.inc"
 
C
      DIMENSION PAR(100)
      REAL PX,PY,PZ,PHIC,PHI1,RMIN,RMAX
      REAL R1,R2,R3,R4,R5,R6,R7,R8,R9
C
C----------------------------------------------------------------------
C
C        Calculates range of each cell
C
      PHI1= (PAR(1)+180)/180*3.14159265359
      NZ= PAR(4)
      NPDV= PAR(3)
      PHIC= (PAR(2)/PAR(3))/180*3.14159265359
C      Rotation around the Z-axis.  Coeffisients of rotation:
      R1= COS(PHI1)
      R2= -SIN(PHI1)
      R3= 0.0
      R4= SIN(PHI1)
      R5= COS(PHI1)
      R6= 0.0
      R7= 0.0
      R8= 0.0
      R9= 1.0
C      Face defined in the yz-plane (x=0)
C----------------------------------------------------------------------
C
      WRITE(BLKSTR,10000)N1
      CALL GJWRIT
C---------------------------------------------------------
C   1. Definition point for the face:
C
      N3= 2
      PX= 0.0
C
      PY= PAR(6)
      PZ= PAR(5)
C---------------------------------------------------------
C
      WRITE(BLKSTR,10100)PX,PY,PZ
      CALL GJWRIT
C---------------------------------------------------------
C  Loops over the other definition points:
C
      DO 10  K=1,NZ
         N3=N3+3
         PY= PAR(N3+2)
         PZ= PAR(N3)
         RMIN= PAR(N3+1)
         RMAX= PAR(N3+2)
         IF (RMIN .GE. RMAX) GOTO 10
C
C
C---------------------------------------------------------
         WRITE(BLKSTR,10100)PX,PY,PZ
         CALL GJWRIT
C---------------------------------------------------------
C
   10 CONTINUE
C
      DO 20  L=2,NZ
         PY= PAR(N3+1)
         PZ= PAR(N3)
C
C---------------------------------------------------------
         WRITE(BLKSTR,10100)PX,PY,PZ
         CALL GJWRIT
C----------------------------------------------------------
         N3=N3-3
   20 CONTINUE
C
C      Geometric Transformation
C*       WRITE SET @302,N1..#301,R1,R2,R3,R4,R5,R6,R7,R8,R9
C
      N1=N1+1
      NG= N1
      WRITE(BLKSTR,10200)N1,R1,R2,R3,R4,R5,R6,R7,R8,R9
      CALL GJWRIT
C-------------------------------------------------
C---First FACE:
C*       WRITE SET @113,F2..#101,!F1,N1
C------------------------------------------------
C
      N1=N1+1
      WRITE(BLKSTR,10300)N1,N1-2,N1-1
      CALL GJWRIT
C---------------------------------------------------------
C      Next step is to obtain the right position of the second face to
C      create a cell:
C      Rotation around the Z-axis.  Coeffisients of rotation:
      R1= COS(PHIC)
      R2= -SIN(PHIC)
      R3= 0.0
      R4= SIN(PHIC)
      R5= COS(PHIC)
      R6= 0.0
      R7= 0.0
      R8= 0.0
      R9= 1.0
C
C      Geometric Transformation
C*       WRITE SET @302,N1..#301,R1,R2,R3,R4,R5,R6,R7,R8,R9
C
C------------------------------------------------------------
      N1=N1+1
      NG= N1
      WRITE(BLKSTR,10200)N1,R1,R2,R3,R4,R5,R6,R7,R8,R9
      CALL GJWRIT
C--------------------------------------------------------------
C      Second Face:
C
C      The first cell of the PGON is defined as a ruled solid between two
C      faces.  The second face is defined as a copy of the first face that is
C      rotated PHIC degrees around the Z-axis.
C
C*       WRITE SET @113,F2..#101,!F1,N1
C------------------------------------------------
C
      N1=N1+1
      WRITE(BLKSTR,10300)N1,N1-2,N1-1
      CALL GJWRIT
C
C-------------------------------------------------
C*       WRITE SET @100,N2..#145,!F1,!F2
C------------------------------------------------
C
      N1=N1+1
      WRITE(BLKSTR,10400)N1,N1-3,N1-1
      CALL GJWRIT
C
C--------------------------------------------------------------
C  The rest of the cells are defined as rotated copies of the first cell:
      N2=N1
      DO 30  K=2, NPDV
C*       WRITE SET @100,N3..#101,!N1-1,!NG
C
C------------------------------------------------
C
         N1=N1+1
         WRITE(BLKSTR,10500)N1,N1-1,NG
         CALL GJWRIT
C
   30 CONTINUE
C------------------------------------------------------------
C      The final shape is a Boolean union of the cells
C*       WRITE SET @100,N4..#100,2,!N1-NPDV....!N1-1
C------------------------------------------------
C
      N1=N1+1
      WRITE(BLKSTR,10700)N1
      CALL GJWRIT
C------------------------------------------------
      DO 40  K=N2, N1-1
C
         WRITE(BLKSTR,10600)K
         CALL GJWRIT
C
   40 CONTINUE
C
10000   FORMAT('@103,',I10,',:5,2#3,3,2')
10100   FORMAT(',',G14.7,',',G14.7,',',G14.7)
10200   FORMAT('@302,',I10,',:5,2#301,',G14.7,',',G14.7,',',G14.7
     +  ,',',G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7)
10300   FORMAT('@113,',I10,',:5,2#101,!',I10,',!',I10)
10400   FORMAT('@100,',I10,',:5,2#145,!',I10,',!',I10)
10500   FORMAT('@100,',I10,',:5,2#101,!',I10,',!',I10)
10600   FORMAT(',!',I10)
10700   FORMAT('@100,',I10,',:5,2#100,2')
C
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDXYZ(IT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draws track IT for which space points have been stored   *
C.    *       in banks JXYZ                                            *
C.    *       If IT=0 draws all tracks                                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Authors : R.Brun, P.Zanarini, S.Giani                    *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      PARAMETER (NTRTYP=8)
      INTEGER COLO(NTRTYP), MODE(NTRTYP)
      DIMENSION U(50),V(50)
C
C             MODE(ITRTYP) and COLO(ITRTYP) are used by ISLN and GDCOL
C             ITRTYP (the track type) correspond to :
C
C             ITRTYP  PARTICLE  |       HIGZ         |
C             ----------------------------------------
C                 1   gammas    | dotted      blue   |
C                 2   electrons | solid       red    |
C                 3   neutral   | dot-dashed  black  |
C                 4   hadrons   | solid       red    |
C                 5   muons     | dashed      green  |
C                 6   geantino  | dot-dashed  black  |
C                 7   cerenkov  | dotted      yellow |
C                 8   ions      | solid       light b|
C
      SAVE MODE,COLO
      DATA MODE /3, 1, 4, 1, 2, 1, 3, 1/
      DATA COLO/4, 2, 1, 2, 3, 1, 6, 5/
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('OFF ',IOFF,4,4)
      CALL UCTOH('LINE',ILINE,4,4)
C
C            Set IOBJ to TRACK
C
      IOBJ=2
C
      IF(JXYZ.LE.0)GO TO 70
      NXYZ=IQ(JXYZ-2)
      N1  =1
      N2  =NXYZ
      IF(IT.NE.0)N1=IT
      IF(IT.NE.0)N2=IT
      IF(N1.LE.0)GO TO 70
      IF(N2.GT.NXYZ)GO TO 70
C
      DO 60 I=N1,N2
         JX=LQ(JXYZ-I)
         IF(JX.LE.0)GO TO 50
         IPOINT=JX+3
   10    IF(IPOINT.GT.JX+IQ(JX-1))GO TO 50
         NPOINT=Q(IPOINT)
         IF(NPOINT.LE.0)GO TO 50
         IPART=Q(IPOINT+1)
         JPA=LQ(JPART-IPART)
         ITRTYP=Q(JPA+6)
C
C             Line style and color code
C
         IMOD=MODE(ITRTYP)
         ICOL=COLO(ITRTYP)
C
         CALL GDCOL(ICOL)
         CALL GDCOL1(ICOL)
C
         IPF=1
   20    IIPL=IPF+49
         IF(IIPL.GT.NPOINT)IIPL=NPOINT
         NP=IIPL-IPF+1
         IF(NP.LE.0)GO TO 40
C
C            Is THRZ option set on (R-Z projection)
C            and the cut lateral (vertical or horizontal) ?
C
         IF (ITHRZ.NE.IOFF) THEN
            IF (ICUT.EQ.1.OR.ICUT.EQ.2) ITR3D=-I
         ELSE
            ITR3D=0
         ENDIF
C
         CALL GDFR3D(Q(IPOINT+3*IPF-1),NP,U,V)
C
         IF (IDVIEW.EQ.0) THEN
C
            IF (ITRKOP.EQ.ILINE.AND.NP.GT.1) THEN
               CALL ISLN(IMOD)
               CALL IGPID(1,'Kine',I,' ')
               CALL IPL(NP,U,V)
            ELSE
               DO 30 IP=1,NP
                  CALL IGPID(1,'Kine',I,' ')
                  CALL IPM(1,U(IP),V(IP))
   30          CONTINUE
            ENDIF
C
         ELSE
C
            CALL IGPID(1,'Kine',I,' ')
            CALL GDRAWP(U,V,NP)
C
         ENDIF
C
         IPF=IIPL
         IF(IPF.LT.NPOINT)GO TO 20
   40    IPOINT=IPOINT+3.*Q(IPOINT)+2.
         GO TO 10
C
   50    CONTINUE
C
   60 CONTINUE
C
      CALL GDCOL(0)
      CALL GDCOL1(IBITS(LINATP,16,4))
C
   70 CONTINUE
C
C             Reset line style
C
      CALL ISLN(1)
C
C             Reset ITR3D to 'standard projection'
C
      ITR3D=0
      IOBJ=0
C
      END

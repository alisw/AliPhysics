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
      SUBROUTINE GDTRAK(CHOPTT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Store space point VECT in banks JXYZ                     *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Authors    R.Brun, F.Carena  *********                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcpush.inc"
#include "geant321/gckine.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
      PARAMETER (NPOBUF=50)
      CHARACTER*(*) CHOPTT
      CHARACTER*10  CHOPT
      LOGICAL OPTDRW, OPTREM, FIRST
      INTEGER COLO(8), OLDCOL, OLDLIN
      DIMENSION U(NPOBUF),V(NPOBUF),MODE(8)
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
C                 7   Cerenkov  | dotted      yellow |
C                 8   Ions      | solid       light b|
C
      SAVE MODE,COLO,FIRST,IOFF,ILINE,OLDCOL,OLDLIN
      DATA MODE /3, 1, 4, 1, 2, 1, 3, 1/
      DATA COLO/4, 2, 1, 2, 3, 1, 6, 5/
      DATA FIRST /.TRUE./
      DATA OLDCOL/ -1 /, OLDLIN / -1 /
*
*     ------------------------------------------------------------------
*
      IF(JXYZ.EQ.0)THEN
*
*  ***         No bank...
*
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
         GOTO 999
      ENDIF
*
      IF(ITRA.GT.IQ(JXYZ-2)) GOTO 40
*
      JX=LQ(JXYZ-ITRA)
      IF(JX.EQ.0) GOTO 40
*
      CHOPT  = CHOPTT
      CALL CLTOU(CHOPT)
      OPTDRW = INDEX(CHOPT,'D') .NE. 0
      OPTREM = INDEX(CHOPT,'R') .NE. 0
      IF(.NOT.(OPTDRW.OR.OPTREM)) THEN
         OPTDRW = .TRUE.
         OPTREM = .TRUE.
      ENDIF
*
      IF (FIRST) THEN
         CALL UCTOH('OFF ',IOFF,4,4)
         CALL UCTOH('LINE',ILINE,4,4)
         FIRST  = .FALSE.
      ENDIF
*
      IF (OPTDRW) THEN
*
*  ***        Draw the current track
*
         IPOINT = JX+Q(JX+2)
         NPOINT = Q(IPOINT)
         IPART  = Q(IPOINT+1)
         JPA    = LQ(JPART-IPART)
         ITRTYP = Q(JPA+6)
*
*  ***        Line style and color code
*
         IMOD   = MODE(ITRTYP)
         ICOL   = COLO(ITRTYP)
         IF(ICOL.NE.OLDCOL) THEN
            CALL GDCOL(ICOL)
            CALL GDCOL1(ICOL)
            OLDCOL=ICOL
         ENDIF
*
         NDONE  = 0
   10    NDO    = MIN(NPOINT-NDONE,NPOBUF)
         IF (NDO.GT.0) THEN
*
*            Is THRZ option set on (R-Z projection)
*            and the cut lateral (vertical or horizontal) ?
*
            IF (ITHRZ.NE.IOFF) THEN
               IF (ICUT.EQ.1.OR.ICUT.EQ.2) ITR3D=-ITRA
            ELSE
               ITR3D=0
            ENDIF
*
            CALL GDFR3D(Q(IPOINT+3*NDONE+2),NDO,U,V)
*
            IF (IDVIEW.EQ.0) THEN
               IF (ITRKOP.EQ.ILINE.AND.NDO.GT.1) THEN
                  IF(OLDLIN.NE.IMOD) THEN
                     CALL ISLN(IMOD)
                     OLDLIN=IMOD
                  ENDIF
                  CALL IPL(NDO,U,V)
               ELSE
                  DO 20 IP=1,NDO
                     CALL IPM(1,U(IP),V(IP))
   20             CONTINUE
               ENDIF
*
            ELSE
               CALL GDRAWP(U,V,NDO)
*
            ENDIF
            NDONE = NDONE + NDO
            GOTO 10
         ENDIF
      ENDIF
      IF (OPTREM) THEN
*
* ***      Remove the current track
*
         IPOINT = JX+INT(Q(JX+2))
         Q(JX+1)= Q(JX+1) + 3.*Q(IPOINT)
         Q(IPOINT) = 0.
         IQ(JX) = -77
      ENDIF
      GOTO 999
*
   40 WRITE (CHMAIL,10100) ITRA
      CALL GMAIL(0,0)
*
10000 FORMAT(' *** GDTRAK: No JXYZ bank')
10100 FORMAT(' *** GDTRAK: No points stored for track ',I10)
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCXYZ
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Online display of the track trajectory                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Authors : R.Brun, P.Zanarini, S.Giani  ****              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gckine.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctrak.inc"
      PARAMETER (NTRTYP=8)
      DIMENSION UU(2),VV(2)
      DIMENSION MODE(NTRTYP),KCOLO(NTRTYP),VOLD(3)
      SAVE MODE,KCOLO,IFIRST,IOFF,ILINE,UU,VV
      DATA MODE /3, 1, 4, 1, 2, 1, 1, 1/
      DATA KCOLO/4, 2, 1, 2, 3, 1, 6, 6/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
*      IF(SLENG.GT.0..AND.STEP.LT.0.1)RETURN
      IF(IFIRST.EQ.0)THEN
         CALL UCTOH('OFF ',IOFF,4,4)
         CALL UCTOH('LINE',ILINE,4,4)
         IFIRST=1
      ENDIF
C
C             Set IOBJ to TRACK
C
      IOBJ=2
C
C             Is THRZ option set on (R-Z projection)
C             and the cut lateral (vertical or horizontal) ?
C
      IF (ITHRZ.NE.IOFF) THEN
         IF (ICUT.EQ.1.OR.ICUT.EQ.2) ITR3D=-ITRA
      ELSE
         ITR3D=0
      ENDIF
      ICOL=KCOLO(ITRTYP)
      KCOL=LOOKTB(ICOL)
C
      IF(NJTMAX.EQ.0)THEN
*
*            Normal sequential tracking
         IF(SLENG.LE.0.)THEN
            CALL GDFR3D(VECT,1,UU(1),VV(1))
            IF (ITRKOP.NE.ILINE) THEN
               CALL IGPID(1,'Kine',ITRA,' ')
               CALL IPM(1,UU,VV)
            ENDIF
         ELSE
            CALL GDFR3D(VECT,1,UU(2),VV(2))
            IF (ITRKOP.EQ.ILINE) THEN
               CALL ISPLCI(KCOL)
               CALL ISLN(MODE(ITRTYP))
               CALL IGPID(1,'Kine',ITRA,' ')
               CALL IPL(2,UU,VV)
            ELSE
               CALL IGPID(1,'Kine',ITRA,' ')
               CALL IPM(1,UU(2),VV(2))
            ENDIF
            UU(1)=UU(2)
            VV(1)=VV(2)
         ENDIF
      ELSE
*
*            Parallel tracking
         CALL GDFR3D(VECT,1,UU(2),VV(2))
         IF (ITRKOP.NE.ILINE) THEN
            CALL IGPID(1,'Kine',ITRA,' ')
            CALL IPM(1,UU(2),VV(2))
         ELSE
            CALL ISPLCI(KCOL)
            CALL ISLN(MODE(ITRTYP))
            VOLD(1)=VECT(1)-STEP*VECT(4)
            VOLD(2)=VECT(2)-STEP*VECT(5)
            VOLD(3)=VECT(3)-STEP*VECT(6)
            CALL GDFR3D(VOLD,1,UU(1),VV(1))
            CALL IGPID(1,'Kine',ITRA,' ')
            CALL IPL(2,UU,VV)
            IF (MODE(ITRTYP).NE.1) CALL ISLN(1)
         ENDIF
      ENDIF
C
C             Reset ITR3D to 'standard projection'
C
      ITR3D=0
      IOBJ=0
C
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTHRZ(XYZ)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine processing R-Z projection of tracks and hits     *
C.    *                                                                *
C.    *       XYZ    = Array of X,Y,Z coordinates (input/output)       *
C.    *                                                                *
C.    *    ==>Called by : GDFR3D                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcnum.inc"
      DIMENSION XYZ(3),ITRSGN(20),ITRFLG(20)
      LOGICAL BTEST
      SAVE MAXITR,EPS,LIDEVT,LITR3D,ITRSGN,ITRFLG,IFACT
      DATA MAXITR/20/,EPS/0.0001/,LIDEVT/0/
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('360 ',I360,4,4)
C
C             If ITHRZ='ON  ' or '180 ' then use positive or negative ra
C             else if ITHRZ='360 ' then use only positive radius
C
      IF (ITHRZ.EQ.I360) THEN
         XYZ(2)=SQRT(XYZ(2)*XYZ(2)+XYZ(3)*XYZ(3))
         GO TO 999
      ENDIF
C
C             IDEVT  = Event number (in common GCFLAG)
C             LIDEVT = Last event number
C
      IF (IDEVT.NE.LIDEVT) THEN
         DO 10 I=1,MAXITR
            ITRSGN(I)=0
            ITRFLG(I)=0
   10    CONTINUE
         LIDEVT=IDEVT
      ENDIF
C
C             NBIT   = Number of bit per machine word (in common GCNUM)
C             MAXITR = A maximum of NBIT*MAXITR significant track signs
C                      can be stored (after that all signs will be posit
C             IWORD  = 1,2,...,MAXITR
C             IOFSET = 0,1,...,NBIT-1
C
      IWORD=(IABS(ITR3D)/NBIT)+1
      IOFSET=MAX(MOD(IABS(ITR3D),NBIT)-1,0)
C
      IF (IWORD.GT.MAXITR) THEN
         IFACT=1
         XYZ(2)=IFACT*SQRT(XYZ(2)*XYZ(2)+XYZ(3)*XYZ(3))
         GO TO 999
      ENDIF
C
C             If being drawing an hit, with sign already
C             set by previous track, then use that sign
C
      IF (IOBJ.EQ.3.AND.BTEST(ITRFLG(IWORD),IOFSET)) THEN
         IF(BTEST(ITRSGN(IWORD),IOFSET)) THEN
            IFACT=1
         ELSE
            IFACT=-1
         ENDIF
         XYZ(2)=IFACT*SQRT(XYZ(2)*XYZ(2)+XYZ(3)*XYZ(3))
         GO TO 999
      ENDIF
C
      IF (ITR3D.LT.0) THEN
C
C             First call to GDFR3D/GDTHRZ : initialize LITR3D and make I
C
         ITR3D=-ITR3D
         IF (.NOT.BTEST(ITRFLG(IWORD),IOFSET)) THEN
            ITRFLG(IWORD) = IBSET(ITRFLG(IWORD),IOFSET)
            LITR3D=0
         ENDIF
      ENDIF
C
      IF (ITR3D.GT.0) THEN
C
C             LITR3D contains the last ITR3D value,
C             or zero if last point was at U-V origin (0,0,...)
C
         IF (ITR3D.EQ.LITR3D) THEN
C
C             Old track : get R sign from ITRSGN into IFACT
C
            IF(BTEST(ITRSGN(IWORD),IOFSET)) THEN
               IFACT=1
            ELSE
               IFACT=-1
            ENDIF
C
         ELSE
C
C             New track : set R sign into ITRSGN and IFACT
C
            IF (ABS(XYZ(2)).GT.EPS) THEN
               IF (XYZ(2)/ABS(XYZ(2)).GT.0.) THEN
                  ITRSGN(IWORD) = IBSET(ITRSGN(IWORD),IOFSET)
                  IFACT=1
               ELSE
                  ITRSGN(IWORD) = IBCLR(ITRSGN(IWORD),IOFSET)
                  IFACT=-1
               ENDIF
            ENDIF
C
            LITR3D=ITR3D
            IF ((ABS(XYZ(1)).LT.EPS).AND.(ABS(XYZ(2)).LT.EPS)) LITR3D=0
C
         ENDIF
C
   20    CONTINUE
         XYZ(2)=IFACT*SQRT(XYZ(2)*XYZ(2)+XYZ(3)*XYZ(3))
C
      ENDIF
  999 RETURN
      END

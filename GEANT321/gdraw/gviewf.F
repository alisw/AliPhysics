*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GVIEWF(U,V,NP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the polyline described by U and V vectors,          *
C.    *       of length NP.                                            *
C.    *                                                                *
C.    *       Depending on IDVIEW it draws on screen (IDVIEW=0)        *
C.    *       or stores in the current view bank (IDVIEW>0).           *
C.    *                                                                *
C.    *       In LINATT (common GCDRAW) there is                       *
C.    *       a bit mask for the line attributes :                     *
C.    *                                                                *
C.    *       Bit  1- 7 = Used by view bank (LENGU)                    *
C.    *       Bit  8-10 = Line width                                   *
C.    *       Bit 11-13 = Line style                                   *
C.    *       Bit 14-16 = Fill area                                    *
C.    *       Bit 17-24 = Line color                                   *
C.    *                                                                *
C.    *    ==>Called by : GDRAWV                                       *
C.    *       Author : P.Zanarini ; S.Giani  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcspee.inc"
      DIMENSION U(*),V(*)
      SAVE LFILOL
      DATA LFILOL/-1/
C.
C.    ------------------------------------------------------------------
C.
C
C             Store on view bank IDVIEW
C
      JV=LQ(JDRAW-IDVIEW)
      IGU=IGU+1
C
   10 IF (IGU.LE.MAXGU) GO TO 20
C
C             Push graphic unit banks
C
      IF(MORGU.EQ.0)MORGU=100
      MORPUS=MAX(MORGU,MAXGU/4)
      JV = LQ(JV-1)
      CALL MZPUSH(IXCONS,JV,0,MORPUS,'I')
      IF(IEOTRI.NE.0)GO TO 40
      JV=LQ(JDRAW-IDVIEW)
      JV = LQ(JV-2)
      CALL MZPUSH(IXCONS,JV,0,MORPUS,'I')
      IF(IEOTRI.NE.0)GO TO 40
      JV=LQ(JDRAW-IDVIEW)
      MAXGU=MAXGU+MORPUS
      GO TO 10
C
   20 IF ((IGS+NP).LE.MAXGS) GO TO 30
C
C             Push graphic segment banks
C
      IF(MORGS.EQ.0)MORGS=100
      MORPUS=MAX(MORGS,MAXGS/4,NP)
      JV = LQ(JV-4)
      CALL MZPUSH(IXCONS,JV,0,MORPUS,'I')
      IF(IEOTRI.NE.0)GO TO 40
      JV=LQ(JDRAW-IDVIEW)
      JV = LQ(JV-5)
      CALL MZPUSH(IXCONS,JV,0,MORPUS,'I')
      IF(IEOTRI.NE.0)GO TO 40
      JV=LQ(JDRAW-IDVIEW)
      MAXGS=MAXGS+MORPUS
*
      GO TO 20
C
   30 CONTINUE
      Q(JV+13)=GTHETA
      Q(JV+14)=GPHI
      Q(JV+15)=GPSI
      Q(JV+16)=GU0
      Q(JV+17)=GV0
      Q(JV+18)=GSCU
      Q(JV+19)=GSCV
      JV1=LQ(JV-1)
      JV2=LQ(JV-2)
      JV4=LQ(JV-4)
      JV5=LQ(JV-5)
*
      CALL UCOPY(U,Q(JV4+IGS+1),NP)
      CALL UCOPY(V,Q(JV5+IGS+1),NP)
C
C             Bit 1- 7 = LENGU
C             Bit 8-24 = Line attribute
C
      ISUM=0
      CALL MVBITS(LINATT,0,24,ISUM,0)
      CALL MVBITS(NP,0,10,ISUM,0)
      Q(JV1+IGU)=ISUM
C
      Q(JV2+IGU)=IGS+1
      IGS=IGS+NP
      GO TO 999
C
   40 WRITE (CHMAIL,10000)
      CALL GMAIL(0,0)
C
10000 FORMAT (' *** GVIEWF ***: Memory overflow in pushing a bank')
  999 END

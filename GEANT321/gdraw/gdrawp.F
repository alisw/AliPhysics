*
* $Id$
*
* $Log$
* Revision 1.2  1996/09/30 13:37:32  ravndal
* Backward compatibility for view banks
*
* Revision 1.1.1.1  1995/10/24 10:20:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRAWP(U,V,NP)
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
C.    *       Author : P.Zanarini ; S.Giani  1992 ********             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcspee.inc"
*
      COMMON/SP3D/ISPFLA
*
      DIMENSION U(*),V(*)
      SAVE LFILOL
      DATA LFILOL/-1/
C.
C.    ------------------------------------------------------------------
C.
      LLEP=ABS(LEP)
      LINFLA=0
      IF (IDVIEW.EQ.0.OR.IDVIEW.EQ.-175) GO TO 40
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
      IF(IEOTRI.NE.0)GO TO 50
      JV=LQ(JDRAW-IDVIEW)
      JV = LQ(JV-2)
      CALL MZPUSH(IXCONS,JV,0,MORPUS,'I')
      IF(IEOTRI.NE.0)GO TO 50
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
      IF(IEOTRI.NE.0)GO TO 50
      JV=LQ(JDRAW-IDVIEW)
      JV = LQ(JV-5)
      CALL MZPUSH(IXCONS,JV,0,MORPUS,'I')
      IF(IEOTRI.NE.0)GO TO 50
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
      IFIL=IBITS(ISUM,13,3)
      IF(IFIL.EQ.0)THEN
         CALL MVBITS(NP,0,7,ISUM,0)
      ELSE
         CALL MVBITS(NP,0,10,ISUM,0)
      ENDIF
      Q(JV1+IGU)=ISUM
C
      Q(JV2+IGU)=IGS+1
      IGS=IGS+NP
      GO TO 999
C
C             Draw vectors on screen
C
   40 CONTINUE
C
C             Extract the new line attributes
C
      LINCOL=IBITS(LINATT,16,8)
      CALL ISFACI(LINCOL)
      LINFIL=IBITS(LINATT,13,3)
      IF(IDVIEW.NE.-175.OR.LINFIL.EQ.0)THEN
         LINWID=IBITS(LINATT,7,3)
         IF(LINWID.GT.1)LINWID=LINWID*2
      ELSE
         LINWID=8-LINFIL
         IF(LINFIL.EQ.1)LINWID=2
         IF(LINWID.GT.1)LINWID=LINWID*2
         IF(ZZFV.GT.1.)LINWID=LINWID*ZZFV
      ENDIF
      LINSTY=IBITS(LINATT,10,3)
      IF(LINSTY.EQ.7)LINSTY=1
      IF(LINFIL.LE.1.OR.IDVIEW.EQ.-175.OR.ISPFLA.EQ.1)
     +CALL ISPLCI(LINCOL)
      WLINW=LINWID
      CALL IGSET('LWID',WLINW)
C
C             If NP=1 draw a marker
C
      IF (NP.EQ.1) THEN
         CALL IPM(1,U,V)
      ELSE
C
C
C             Fill area
C
*SG
         IF(IDVIEW.EQ.-175)THEN
            IF(LINFIL.GT.0.AND.NP.GT.2.AND.LINSTY.NE.6)THEN
               CALL ISFAIS(1)
               CALL IFA(NP,U,V)
            ENDIF
         ENDIF
C
C             If NP>1 draw a line with a given style
C             and draw black edges both for HIDE OFF
C             and SHAD options in case of FILL
C
         CALL UCTOH('ON  ',IFLH,4,4)
         IF(IHIDEN.NE.IFLH.AND.LINFIL.GT.0)THEN
           CALL ISPLCI(1)
         ENDIF
         IF(LINSTY.EQ.6.AND.LINFIL.NE.0)THEN
            LINSTY=1
            LINFLA=1
            CALL ISPLCI(1)
            CALL IGSET('LWID',3.)
            IF(LINWID.GE.12)CALL IGSET('LWID',6.)
            IF(LINWID.LE.4)CALL IGSET('LWID',1.)
         ENDIF
         IF(LLEP.LE.10.OR.LINFIL.EQ.0.OR.LINFLA.NE.1)THEN
            IF (LINSTY.EQ.1) THEN
***            call write_dxf_pline(np,u,v,lincol,linwid,1)
               CALL IPL(NP,U,V)
C
            ELSE IF (LINSTY.GT.1.AND.LINSTY.LE.4) THEN
               CALL ISLN(LINSTY)
               CALL IPL(NP,U,V)
               CALL ISLN(1)
C
            ENDIF
         ENDIF
C
      ENDIF
C
      GO TO 999
C
   50 WRITE (CHMAIL,10000)
      CALL GMAIL(0,0)
C
10000 FORMAT (' *** GDRAWP ***: Memory overflow in pushing a bank')
  999 END


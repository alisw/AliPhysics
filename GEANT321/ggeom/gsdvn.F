*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSDVN(KNAME,MOTHER,NDIV,IAXIS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine divides MOTHER into NDIV divisions called NAME   *
C.    *       along axis number IAXIS.                                 *
C.    *          JVO=Pointer to MOTHER volume                          *
C.    *          JDIV=LQ(JVO-1)                                        *
C.    *                                                                *
C.    *            Q(JDIV+1)=IAXIS                                     *
C.    *            Q(JDIV+2)=Volume number.                            *
C.    *            Q(JDIV+3)=NDIV                                      *
C.    *            Q(JDIV+4)=Lowest coord of slices.                   *
C.    *            Q(JDIV+5)=Step size in coordinates.                 *
C.    *                                                                *
C.    *    ==>Called by :  <USER>, GEDITV                              *
C.    *         Authors R.Brun,   A.McPherson  *********               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcshno.inc"
      CHARACTER*4 KNAME,MOTHER
      DIMENSION PAR(100),PARM(100),ATT(20)
      SAVE ATT
      DATA ATT /1.,1.,1.,1.,1.,15*0./
C.
C.    ------------------------------------------------------------------
C.
C              Check if volume master bank exists.
C
      CALL UCTOH(KNAME,NAME,4,4)
      IF(JVOLUM.GT.0)GO TO 10
      WRITE(CHMAIL,1000)
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if MOTHER volume exists.
C
  10  CALL GLOOK(MOTHER,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.GT.0)GO TO 20
      WRITE(CHMAIL,1100)MOTHER
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if NAME volume exists.
C
  20  CALL GLOOK(KNAME,IQ(JVOLUM+1),NVOLUM,IN)
      IF(IN.LE.0)GO TO 50
      WRITE(CHMAIL,2000)NAME
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if MOTHER is not divided.
C
  50  JVO=LQ(JVOLUM-IVO)
      NIN=Q(JVO+3)
      IF(NIN.EQ.0)GO TO 60
      WRITE(CHMAIL,4000)MOTHER
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check validity of axis value.
C
  60  IF(IAXIS.GT.0.AND.IAXIS.LT.4)GO TO 70
      WRITE(CHMAIL,5000)IAXIS
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check validity of NDIV
C
  70  IF(NDIV.GT.0)GO TO 80
      WRITE(CHMAIL,6000)NDIV
      CALL GMAIL(0,0)
      GO TO 99
C
C               Create bank to store division parameters.
C
  80  CALL MZBOOK(IXCONS,JDIV,JVO,-1,'VODI',0,0,6,3,0)
      IF(IEOTRI.NE.0)GO TO 95
      IQ(JDIV-5)=IVO
C
C               Now store parameters into bank area.
C
  90  Q(JDIV+1)=IAXIS
      Q(JDIV+2)=NVOLUM+1
      Q(JDIV+3)=NDIV
      Q(JVO+3)=-1
      IVOM= IVO
      NWM = IQ(JVO-1)
      NW  = NWM
      ISH = Q(JVO+2)
C
C               Bit to allow division of objects defined
C               by GSPOSP.
C
      C0=0.0
      STEP=0.0
      NPAR=Q(JVO+5)
      NATT=Q(JVO+6)
      CALL UCOPY(Q(JVO+NPAR+7),ATT,NATT)
      IF(NPAR.LE.0) GO TO 230
C
      CALL GFIPAR(JVO,0,0,NPAR,NATT,PAR,ATT)
      CALL UCOPY(PAR,PARM,NPAR)
C
C              Find and store start and step.
C
      IF(ISH.NE.1) GO TO 100
C
C               Box.
C
      STEP=-1.0
      PAR(IAXIS)=-1.0
      IF(PARM(IAXIS).LT.0.0) GO TO 230
      C0  =-PARM(IAXIS)
      STEP=PARM(IAXIS)*2.0/NDIV
      PAR(IAXIS)=STEP/2.
      GO TO 230
C
  100 CONTINUE
      IF(ISH.NE.2) GO TO 110
C
C              Trapezoid with only X thickness varying with Z.
C
      IF(IAXIS.EQ.1) GO TO 900
      PAR(1)=-1.
      PAR(2)=-1.
      STEP=-1.0
      PAR(IAXIS+1)=-1.0
      IF(PARM(IAXIS+1).LT.0.0) GO TO 230
      C0  =-PARM(IAXIS+1)
      STEP=PARM(IAXIS+1)*2.0/NDIV
      PAR(IAXIS+1)=STEP/2.
      GO TO 230
C
  110 CONTINUE
      IF(ISH.NE.3) GO TO 120
C
C              Trapezoid with both X and Y thicknesses varying with
C              Z
C
      IF(IAXIS.NE.3) GO TO 900
      PAR(1)=-1.
      PAR(2)=-1.
      PAR(3)=-1.
      PAR(4)=-1.
      STEP=-1.0
      PAR(5)=-1.0
      IF(PARM(5).LT.0.0) GO TO 230
      C0  =-PARM(5)
      STEP=PARM(5)*2.0/NDIV
      PAR(5)=STEP/2.
      GO TO 230
C
  120 CONTINUE
      IF(ISH.NE.4) GO TO 125
      IF(IAXIS.NE.3) GO TO 126
      PAR(1)=-1.
      PAR(4)=-1.
      PAR(5)=-1.
      PAR(6)=-1.
      PAR(8)=-1.
      PAR(9)=-1.
      PAR(10)=-1.
      STEP=-1.0
      IF(PARM(1).LT.0.0) GO TO 230
      C0=-PARM(1)
      STEP=PARM(1)*2.0/NDIV
      PAR(1)=STEP*0.5
C
      GO TO 230
C
  126 IF(IAXIS.NE.2) GO TO 900
      IF(MOD(PARM(3),180.).EQ.0.) GO TO 127
      WRITE(CHMAIL,10100)
10100 FORMAT(' Division of TRAP ',A4,
     +    ' along Y only possible when PHI=0,180')
      CALL GMAIL(0,0)
      GOTO 99
  127 IF(PARM(4).EQ.PARM(8))  GO TO 128
      WRITE(CHMAIL,10200)
10200 FORMAT(' Division of TRAP ',A4,
     +    ' along Y only possible when H1=H2')
      CALL GMAIL(0,0)
      GOTO 99
  128 CONTINUE
      STEP = -1.
      IF(PARM(4).LT.0.0) GO TO 230
      C0=-PARM(4)
      STEPH = PARM(4)/NDIV
      PAR(4) = STEPH
      PAR(5) = -1.
      PAR(6) = -1.
      PAR(8) = STEPH
      PAR(9) = -1.
      PAR(10) = -1.
      STEP = 2.*STEPH
C
      GO TO 230
C
  125 CONTINUE
      IF(ISH.NE.5.AND.ISH.NE.6.AND.ISH.NE.NSCTUB) GO TO 160
C
C              Tube, tube segment or cut tube.
C
      IF(IAXIS.NE.3) GO TO 130
      STEP=-1.0
      PAR(3)=-1.0
      IF(PARM(3).LT.0.0) GO TO 230
      C0  =-PARM(3)
      STEP=PARM(3)*2.0/NDIV
      PAR(3)=STEP/2.
      GO TO 230
C
  130 CONTINUE
      IF(IAXIS.NE.1) GO TO 140
      PAR(1)=-1.
      PAR(2)=-1.
      STEP=-1.0
      IF(PARM(1).LT.0.0) GO TO 230
      C0  =PARM(1)
      IF(PARM(2).LT.0.0) GO TO 230
      STEP=(PARM(2)-PARM(1))/NDIV
      GO TO 230
C
  140 CONTINUE
      IF(ISH.EQ.6) GO TO 150
      NW=NW+2
      ISH=6
      C0  =0.0
      STEP=360.0/NDIV
      NPAR=5
      PAR(4)=-STEP/2.
      PAR(5)=STEP/2.
      GO TO 230
C
  150 CONTINUE
      DP=PAR(5)-PAR(4)
      IF(DP.LT.0.0) DP=DP+360.0
      C0  =PAR(4)
      STEP=DP/NDIV
      PAR(4)=-STEP/2.
      PAR(5)=STEP/2.
      GO TO 230
C
  160 CONTINUE
C
      IF(ISH.NE.7.AND.ISH.NE.8) GO TO 190
      IF(IAXIS.EQ.1) GO TO 165
      IF(IAXIS.NE.3) GO TO 170
C
      STEP=-1.0
      PAR(1)=-1.0
      IF(PARM(1).LT.0.0) GO TO 165
      C0=-PARM(1)
      STEP=PARM(1)*2.0/NDIV
      PAR(1)=STEP*0.5
C
  165 CONTINUE
C
      PAR(2)=-1.0
      PAR(3)=-1.0
      PAR(4)=-1.0
      PAR(5)=-1.0
      GO TO 230
C
  170 CONTINUE
      IF(IAXIS.EQ.1) GO TO 230
C
      IF(ISH.EQ.8) GO TO 180
      NW=NW+2
      ISH=8
      C0  =0.0
      STEP=360.0/NDIV
      NPAR=7
      PAR(6)=-STEP/2.
      PAR(7)=STEP/2.
      GO TO 230
C
  180 CONTINUE
      DP=PAR(7)-PAR(6)
      IF(DP.LT.0.0) DP=DP+360.0
      C0  =PAR(6)
      STEP=DP/NDIV
      PAR(6)=-STEP/2.
      PAR(7)=STEP/2.
      GO TO 230
C
  190 CONTINUE
      IF(ISH.NE.9) GO TO 200
      IF(IAXIS.NE.1) GO TO 195
      PAR(1)=-1.0
      PAR(2)=-1.0
      C0    = PARM(1)
      STEP = (PARM(2)-PARM(1))/NDIV
  195 CONTINUE
      IF(IAXIS.NE.2) GO TO 196
      WRITE(CHMAIL,8102)
      CALL GMAIL(0,0)
      GOTO 99
C
  196 CONTINUE
      IF(IAXIS.NE.3) GO TO 230
      C0=PAR(5)
      DP=PAR(6)-PAR(5)
      IF(DP.LE.0.0) DP=DP+360.0
      STEP=DP/NDIV
      PAR(3)=-1.
      PAR(4)=-1.
      PAR(5)=-0.5*STEP
      PAR(6)=0.5*STEP
      GO TO 230
C
  200 CONTINUE
C
      IF(ISH.NE.10) GO TO 210
C
C              Parallelipiped.
C
      C0  =-PAR(IAXIS)
      STEP=-2.0*C0/NDIV
      PAR(IAXIS)=STEP/2.
      GO TO 230
C
  210 CONTINUE
      IF(ISH.GT.12) GO TO 900
      IF(IAXIS.EQ.1) GO TO 230
      IF(IAXIS.EQ.2) GO TO 220
C
      IPNZ=4
      IF(ISH.EQ.12) IPNZ=3
      IF(PAR(IPNZ).NE.2) GO TO 910
C
      ZH=PAR(IPNZ+4)
      ZL=PAR(IPNZ+1)
      STEP=(ZH-ZL)/NDIV
      C0=ZL
      PAR(IPNZ+4)=STEP*0.5
      PAR(IPNZ+1)=-PAR(IPNZ+4)
      PAR(IPNZ+2)=-1.
      PAR(IPNZ+3)=-1.
      PAR(IPNZ+5)=-1.
      PAR(IPNZ+6)=-1.
C
      GO TO 230
  220 CONTINUE
C
      NDV=NDIV
      IF(ISH.EQ.11) NDV=PAR(3)
      Q(JDIV+3)=NDV
      C0=PAR(1)
      STEP=PAR(2)/NDV
      PAR(1)=-STEP*0.5
      PAR(2)=STEP
      IF(ISH.EQ.11)PAR(3)=1.
C
  230 CONTINUE
C
C                Now create the volume for division.
C
      Q(JDIV+4)=C0
      Q(JDIV+5)=STEP
      NVOLUM=NVOLUM+1
      NVOL  =IQ(JVOLUM-2)
      IF(NVOLUM.GT.NVOL)CALL MZPUSH(IXCONS,JVOLUM,50,50,'I')
      CALL MZBOOK(IXCONS,JVO,JVOLUM,-NVOLUM,'VOL1',50,50,NW,3,0)
      IF(IEOTRI.NE.0)GO TO 95
      IQ(JVOLUM+NVOLUM)=NAME
C
C              Copy parameters in data area.
C
      JVOM=LQ(JVOLUM-IVOM)
      CALL UCOPY(Q(JVOM+1),Q(JVO+1),NWM)
      IF(NPAR.GT.0) CALL GSIPAR(JVO,0,NPAR,NATT,PAR,ATT)
      Q(JVO+2)=ISH
      Q(JVO+3)=0.
      GO TO 99
C
  900 CONTINUE
C
C        Divide action not supported.
C
      WRITE(CHMAIL,8000)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,8001) ISH,IAXIS
      CALL GMAIL(0,0)
C
      GO TO 99
C
  910 CONTINUE
C
C          Trying to divide multi Z sector shape along Z.
C
      WRITE(CHMAIL,8100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,8101) ISH,IAXIS,IPNZ,INT(PAR(IPNZ))
      CALL GMAIL(0,0)
C
      GO TO 99
C
C              Not enough space.
C
  95  WRITE(CHMAIL,7000)NAME,MOTHER
      CALL GMAIL(0,0)
C
  99  CONTINUE
 1000 FORMAT(' ***** GSDVN CALLED AND NO VOLUMES DEFINED *****')
 1100 FORMAT(' ***** GSDVN MOTHER VOLUME ',A4,' DOES NOT EXIST *****')
 2000 FORMAT(' ***** GSDVN VOLUME ',A4,' ALREADY EXISTS *****')
 3000 FORMAT(' ***** GSDVN ROTATION MATRIX',I5,' DOES NOT EXISTS *****')
 4000 FORMAT(' ***** GSDVN MOTHER ',A4,' ALREADY DIVIDED *****')
 5000 FORMAT(' ***** GSDVN BAD AXIS VALUE ',I5,' *****')
 6000 FORMAT(' ***** GSDVN BAD NUMBER OF DIVISIONS ',I5,' *****')
 7000 FORMAT(' ***** GSDVN NOT ENOUGH SPACE TO STORE DIVISIONS ',
     +       ' IN ',A4,' *****')
 8000 FORMAT(' DIVIDE ACTION REQUESTED NOT SUPPORTED AT PRESENT.')
 8001 FORMAT(' ISH =',I5,' IAXIS =',I5)
 8100 FORMAT(' ATTEMPT TO DIVIDE MULTI Z SECTOR SHAPE ALONG Z.')
 8101 FORMAT(' ISH =',I5,' IAXIS =',I5,' NZ (THE',I3,'TH PAR) IS',I5)
 8102 FORMAT(' DIVISION OF A SPHERE ALONG AXIS 2 NOT SUPPORTED')
      END

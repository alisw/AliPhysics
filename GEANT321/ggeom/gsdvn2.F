*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSDVN2(KNAME,MOTHER,NDIV,IAXIS,C0I,NUMED)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       DIVIDES MOTHER INTO NDIV DIVISIONS CALLED NAME           *
C.    *       ALONG AXIS IAXIS STARTING AT COORDINATE VALUE C0.        *
C.    *       THE NEW VOLUME CREATED WILL BE MEDIUM NUMBER NUMED.      *
C.    *                                                                *
C.    *          JVO=POINTER TO MOTHER VOLUME                          *
C.    *          JDIV=LQ(JVO-1)                                        *
C.    *                                                                *
C.    *            Q(JDIV+1)=IAXIS                                     *
C.    *            Q(JDIV+2)=VOLUME NUMBER                             *
C.    *            Q(JDIV+3)=NDIV                                      *
C.    *            Q(JDIV+4)=C0                                        *
C.    *            Q(JDIV+5)=STEP SIZE IN COORDINATES.                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GSDVX                                *
C.    *         Authors F.Bruyant,   A.McPherson  *********            *
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
      DIMENSION PAR(50),ATT(20)
      SAVE ATT
      DATA ATT /1.,1.,1.,1.,1.,15*0./
C.
C.    ------------------------------------------------------------------
C.
C              CHECK IF VOLUME MASTER BANK EXISTS
C
      CALL UCTOH(KNAME,NAME,4,4)
      IF(JVOLUM.GT.0)GO TO 10
      WRITE(CHMAIL,1000)
      CALL GMAIL(0,0)
      GO TO 99
C
C              CHECK IF MOTHER VOLUME EXISTS
C
  10  CALL GLOOK(MOTHER,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.GT.0)GO TO 20
      WRITE(CHMAIL,2000)MOTHER
      CALL GMAIL(0,0)
      GO TO 99
C
C              CHECK IF NAME VOLUME EXISTS
C
  20  CALL GLOOK(KNAME,IQ(JVOLUM+1),NVOLUM,IN)
      IF(IN.LE.0)GO TO 50
      WRITE(CHMAIL,2000)NAME
      CALL GMAIL(0,0)
      GO TO 99
C
C              CHECK IF MOTHER IS NOT DIVIDED
C
  50  JVO=LQ(JVOLUM-IVO)
      NIN=Q(JVO+3)
      IF(NIN.EQ.0)GO TO 60
      WRITE(CHMAIL,4000)MOTHER
      CALL GMAIL(0,0)
      GO TO 99
C
C              CHECK VALIDITY OF AXIS VALUE
C
  60  IF(IAXIS.GT.0.AND.IAXIS.LT.4)GO TO 70
      WRITE(CHMAIL,5000)IAXIS
      CALL GMAIL(0,0)
      GO TO 99
C
C              CHECK VALIDITY OF NDIV
C
  70  IF(NDIV.GT.0)GO TO 80
      WRITE(CHMAIL,6000)NDIV
      CALL GMAIL(0,0)
      GO TO 99
C
C               CREATE BANK TO STORE DIVISION PARAMETERS
C
  80  CALL MZBOOK(IXCONS,JDIV,JVO,-1,'VODI',0,0,6,3,0)
      IF(IEOTRI.NE.0)GO TO 95
      IQ(JDIV-5)=IVO
C
C               NOW STORE PARAMETERS INTO BANK AREA
C
  90  Q(JDIV+1)=IAXIS
      Q(JDIV+2)=NVOLUM+1
      Q(JDIV+3)=NDIV
      Q(JVO+3)=-1
      IVOM= IVO
      NWM = IQ(JVO-1)
      NW  = NWM
      CALL GFIPAR(JVO,0,0,NPAR,NATT,PAR,ATT)
C
C              CHECK START AND FIND AND STORE STEP.
C
      ISH=Q(JVO+2)
      C0 = C0I
      IF(ISH.NE.1) GO TO 100
C
C               BOX
C
      IF(PAR(IAXIS).LE.0.0) GO TO 920
      IF(ABS(C0).GT.PAR(IAXIS)) GO TO 910
      STEP=(PAR(IAXIS)-C0)/NDIV
      PAR(IAXIS)=STEP/2.
      GO TO 210
C
  100 CONTINUE
      IF(ISH.NE.2) GO TO 110
C
C              TRAPEZOID WITH ONLY X THICKNESS VARYING WITH Z.
C
      IF(IAXIS.EQ.1) GO TO 900
      IF(PAR(IAXIS+1).LE.0.0) GO TO 920
      IF(ABS(C0).GT.PAR(IAXIS+1)) GO TO 910
      STEP=(PAR(IAXIS+1)-C0)/NDIV
      PAR(1)=-1.
      PAR(2)=-1.
      PAR(IAXIS+1)=STEP/2.
      GO TO 210
C
  110 CONTINUE
      IF(ISH.NE.3) GO TO 120
C
C              TRAPEZOID WITH BOTH X AND Y THICKNESSES VARYING WITH
C              Z
C
      IF(IAXIS.NE.3) GO TO 900
      IF(PAR(5).LE.0.0) GO TO 920
      IF(ABS(C0).GT.PAR(5)) GO TO 910
      STEP=(PAR(5)-C0)/NDIV
      PAR(1)=-1.
      PAR(2)=-1.
      PAR(3)=-1.
      PAR(4)=-1.
      PAR(5)=STEP/2.
      GO TO 210
C
  120 CONTINUE
      IF(ISH.NE.4) GO TO 125
      IF(IAXIS.NE.3) GO TO 126
      IF(PAR(1).LE.0.0) GO TO 920
      IF(ABS(C0).GT.PAR(1)) GO TO 910
      STEP=(PAR(1)-C0)/NDIV
      PAR(1)=STEP*0.5
      PAR(4)=-1.0
      PAR(5)=-1.0
      PAR(6)=-1.0
      PAR(8)=-1.0
      PAR(9)=-1.0
      PAR(10)=-1.0
      GO TO 210
C
  126 IF(IAXIS.NE.2) GO TO 900
      IF(MOD(PAR(3),180.).EQ.0.) GO TO 127
      WRITE(CHMAIL,10100)
10100 FORMAT(' Division of TRAP ',A4,
     +    ' along Y only possible when PHI=0,180')
      CALL GMAIL(0,0)
      GOTO 99
  127 IF(PAR(4).EQ.PAR(8))  GO TO 128
      WRITE(CHMAIL,10200)
10200 FORMAT(' Division of TRAP ',A4,
     +    ' along Y only possible when H1=H2')
      CALL GMAIL(0,0)
      GOTO 99
  128 CONTINUE
      IF(PAR(4).LE.0.) GO TO 920
      IF(ABS(C0).GT.PAR(4)) GO TO 910
      STEP = (PAR(4)-C0)/NDIV
      PAR(4) = 0.5*STEP
      PAR(5) = -1.
      PAR(6) = -1.
      PAR(8) = 0.5*STEP
      PAR(9) = -1.
      PAR(10) = -1.
C
      GO TO 210
C
  125 CONTINUE
      IF(ISH.NE.5.AND.ISH.NE.6.AND.ISH.NE.NSCTUB) GO TO 160
C
C              Tube, tube segment or cut tube.
C
      IF(IAXIS.NE.3) GO TO 130
      IF(PAR(3).LE.0.0) GO TO 920
      IF(ABS(C0).GT.PAR(3)) GO TO 910
      STEP=(PAR(3)-C0)/NDIV
      PAR(3)=STEP/2.
      GO TO 210
C
  130 CONTINUE
      IF(IAXIS.NE.1) GO TO 140
      IF(PAR(1).LE.0.0.OR.PAR(2).LE.0.0) GO TO 920
      IF(C0.LT.PAR(1).OR.C0.GT.PAR(2)) GO TO 910
      STEP=(PAR(2)-C0)/NDIV
      PAR(1)=-1.
      PAR(2)=-1.
      GO TO 210
C
  140 CONTINUE
      IF(ISH.EQ.6) GO TO 150
      NW=NW+2
      ISH=6
      STEP=360.0/NDIV
      NPAR=5
      PAR(4)=-STEP/2.
      PAR(5)=STEP/2.
      GO TO 210
C
  150 CONTINUE
      DP=PAR(5)-PAR(4)
      DC0P = C0-PAR(4)
      SG = SIGN(1.0,DC0P)
      DC0P = MOD( ABS(DC0P), 360.0)
      IF(SG.LE.0.0) DC0P = 360.0-DC0P
      C0 = PAR(4)+DC0P
      IF(C0-PAR(4).LT.0.0) C0=C0+360.0
      IF(C0-PAR(4).GT.DP) GO TO 910
      DP=PAR(5)-C0
      IF(DP.LT.0.0) DP=DP+360
      STEP=DP/NDIV
      PAR(4)=-STEP/2.
      PAR(5)=STEP/2.
      GO TO 210
C
  160 CONTINUE
C
      IF(ISH.NE.7.AND.ISH.NE.8) GO TO 190
      IF(IAXIS.EQ.1) GO TO 165
      IF(IAXIS.NE.3) GO TO 170
C
      IF(PAR(1).LE.0.0) GO TO 920
      IF(ABS(C0).GT.PAR(1)) GO TO 910
      STEP=(PAR(1)-C0)/NDIV
      PAR(1)=STEP*0.5
C
  165 CONTINUE
C
      PAR(2)=-1.0
      PAR(3)=-1.0
      PAR(4)=-1.0
      PAR(5)=-1.0
      GO TO 210
C
  170 CONTINUE
      IF(IAXIS.EQ.1) GO TO 210
C
      IF(ISH.EQ.8) GO TO 180
      NW=NW+2
      ISH=8
      STEP=360.0/NDIV
      NPAR=7
      PAR(6)=-STEP/2.
      PAR(7)=STEP/2.
      GO TO 210
C
  180 CONTINUE
      DP=PAR(7)-PAR(6)
      DC0P = C0-PAR(6)
      SG = SIGN(1.0,DC0P)
      DC0P = MOD( ABS(DC0P), 360.0)
      IF(SG.LE.0.0) DC0P = 360.0-DC0P
      C0 = PAR(6)+DC0P
      IF(C0-PAR(6).LT.0.0) C0=C0+360.0
      IF(C0-PAR(6).GT.DP) GO TO 910
      DP=PAR(7)-C0
      IF(DP.LT.0.0) DP=DP+360.0
      STEP=DP/NDIV
      PAR(6)=-STEP/2.
      PAR(7)=STEP/2.
      GO TO 210
C
  190 CONTINUE
      IF(ISH.NE.9) GO TO 200
      IF(IAXIS.NE.1) GO TO 195
      PAR(1)=-1.0
      PAR(2)=-1.0
      IF(C0.LT.PAR(1).OR.C0.GT.PAR(2)) THEN
         GOTO 910
      ENDIF
      STEP = (PAR(2)-C0)/NDIV
      GOTO 210
  195 CONTINUE
      IF(IAXIS.NE.2) GO TO 196
      WRITE(CHMAIL,8102)
      CALL GMAIL(0,0)
      GOTO 99
C
  196 CONTINUE
      IF(IAXIS.NE.3) GO TO 210
      ANGMIN = MOD(PAR(5),360.)
      IF(ANGMIN.LT.0.) ANGMIN=ANGMIN+360.
      ANGMAX = MOD(PAR(6),360.)
      IF(ANGMAX.LE.ANGMIN) ANGMAX=ANGMAX+360.
      C0 = MOD(C0,360.)
      IF(C0.LT.0.0) C0=C0+360.0
      IF(C0.GT.ANGMAX.OR.C0.LT.ANGMIN) GO TO 910
      STEP=(ANGMAX-C0)/NDIV
      PAR(3)=-1.
      PAR(4)=-1.
      PAR(5)=-0.5*STEP
      PAR(6)=0.5*STEP
      GO TO 210
C
  200 CONTINUE
C
      GO TO 900
C
C                NOW CREATE THE VOLUME FOR DIVISION
C
 210  Q(JDIV+4)=C0
      Q(JDIV+5)=STEP
      NVOLUM=NVOLUM+1
      NVOL  =IQ(JVOLUM-2)
      IF(NVOLUM.GT.NVOL)CALL MZPUSH(IXCONS,JVOLUM,50,50,'I')
      CALL MZBOOK(IXCONS,JVO,JVOLUM,-NVOLUM,'VOL1',50,50,NW,3,0)
      IF(IEOTRI.NE.0)GO TO 95
      IQ(JVOLUM+NVOLUM)=NAME
C
C              COPY PARAMETERS IN DATA AREA
C
      JVOM=LQ(JVOLUM-IVOM)
      CALL UCOPY(Q(JVOM+1),Q(JVO+1),NWM)
      CALL GSIPAR(JVO,0,NPAR,NATT,PAR,ATT)
      Q(JVO+2)=ISH
      Q(JVO+3)=0.
      Q(JVO+4)=NUMED
      GO TO 99
C
  900 CONTINUE
C
C        DIVIDE ACTION NOT SUPPORTED.
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
C              C0 START OF DIVISION OUT OF OBJECT.
C
      WRITE(CHMAIL,9000) C0
      CALL GMAIL(0,0)
C
      GO TO 99
  920 CONTINUE
C
C              +VE DEFINITE PARAMETER IN DIMENSION OF C0 SET -VE OR 0.
C
      WRITE(CHMAIL,9010)
      CALL GMAIL(0,0)
C
      GO TO 99
C
C              NOT ENOUGH SPACE
C
  95  WRITE(CHMAIL,7000)NAME,MOTHER
      CALL GMAIL(0,0)
C
  99  CONTINUE
 1000 FORMAT(' ***** GSDVN2 CALLED AND NO VOLUMES DEFINED *****')
 2000 FORMAT(' ***** GSDVN2 VOLUME ',A4,' ALREADY EXISTS *****')
 3000 FORMAT(' ***** GSDVN2 ROTATION MATRIX',I5,' DOES NOT EXIST *****')
 4000 FORMAT(' ***** GSDVN2 MOTHER ',A4,' ALREADY DIVIDED *****')
 5000 FORMAT(' ***** GSDVN2 BAD AXIS VALUE ',I5,' *****')
 6000 FORMAT(' ***** GSDVN2 BAD NUMBER OF DIVISIONS ',I5,' *****')
 7000 FORMAT(' ***** GSDVN2 NOT ENOUGH SPACE TO STORE DIVISIONS ',
     +       ' IN ',A4,' *****')
 8000 FORMAT(' DIVIDE ACTION WITH C0 REQUESTED NOT SUPPORTED',
     +' AT PRESENT.')
 8001 FORMAT(' ISH =',I5,' IAXIS =',I5)
 8102 FORMAT(' DIVISION OF A SPHERE ALONG AXIS 2 NOT SUPPORTED')
 9000 FORMAT(' ***** GSDVN2 C0',E15.5,' OUT OF OBJECT *****')
 9010 FORMAT(' ***** GSDVN2 C0 WITH -VE DIMENSION IN MOTHER *****')
      END

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
      SUBROUTINE GSATT(NAME,IAT,IVAL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        Set attribute IAT to the value IVAL                     *
C.    *        for the volume NAME                                     *
C.    *                                                                *
C.    *    1:  'WORK'   0=volume inactive   1=active volume            *
C.    *    2:  'SEEN'   0=unseen 1=seen -1,-2=tree unseen -3=tree opt. *
C.    *    3:  'LSTY'   line style 1,2,3,...                           *
C.    *    4:  'LWID'   line width 1,2,3,...                           *
C.    *    5:  'COLO'   color code 1,2,3,...                           *
C.    *    6:  'FILL'   fill area  0,1,2,...                           *
C.    *    7:  'SET '   set number associated to this volume           *
C.    *    8:  'DET '   detector number associated to this volume      *
C.    *    9:  'DTYP'   detector type (1,2)                            *
C.    *   10:  'NODE'   <>0=a node is created for 3D (PIONS,GMR,etc.)  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDPRTR,GDSPEC,GDTREE,GSDET, *
C.    *                   GSDETV,GSDVN,GSDVN2,GSDVT,GSDVT2,GSVOLU      *
C.    *         Author  R.Brun  S.Giani *********                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      COMMON/DEFPAR/IATDEF,IVADEF
*
      CHARACTER*4 KATT(10),NAME,IAT,IATDEF
      DIMENSION LATT(10)
      SAVE KATT,IFIRST,LATT
      DATA KATT/'WORK','SEEN','LSTY','LWID','COLO','FILL','SET ',
     +'DET ','DTYP','NODE'/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
C              CHECK IF OPTION EXISTS
*
      IF(IAT.EQ.'DEFA')IAT=IATDEF
      IF(IVAL.EQ.10000)IVAL=IVADEF
C
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH(KATT,LATT,4,40)
      ENDIF
C
      IPN=5
      IPS=7
C
C              CHECK IF OPTION EXISTS
C
      CALL GLOOK(IAT,LATT,10,IOP)
      IF(IOP.GT.0)GO TO 10
      WRITE(CHMAIL,10100)IAT
      CALL GMAIL(0,0)
      GO TO 999
C              CHECK IF VOLUME EXISTS
C
   10 IF(JVOLUM.GT.0)GO TO 20
      WRITE(CHMAIL,10000)NAME
      CALL GMAIL(0,0)
      GO TO 999
C
   20 IF(NAME(1:1).EQ.'*')GO TO 40
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVOLU)
      IF(IVOLU.GT.0)GO TO 30
      WRITE(CHMAIL,10000)NAME
      CALL GMAIL(0,0)
      GO TO 999
C
   30 JVO=LQ(JVOLUM-IVOLU)
      NPAR=Q(JVO+IPN)
      Q(JVO+NPAR+IPS-1+IOP)=IVAL
      GO TO 999
C
C             SET OPTION VALUE FOR ALL DEFINED VOLUMES
C
   40 CONTINUE
      DO 50 I=1,NVOLUM
         JVO=LQ(JVOLUM-I)
         NPAR=Q(JVO+IPN)
         Q(JVO+NPAR+IPS-1+IOP)=IVAL
   50 CONTINUE
C
10000 FORMAT(' GSATT VOLUME ',A4,' DOES NOT EXIST ')
10100 FORMAT(' GSATT OPTION ',A4,' NOT IMPLEMENTED ')
  999 END

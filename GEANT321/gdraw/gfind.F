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
      SUBROUTINE GFIND(NAME,IAT,ISVALS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        For attribute IAT find value ISVALS                     *
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
C.    *          Author  S.Giani  *********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 KATT(10),NAME,IAT
      DIMENSION LATT(10)
      SAVE KATT,IFIRST,LATT
      DATA KATT/'WORK','SEEN','LSTY','LWID','COLO','FILL','SET ',
     +'DET ','DTYP','NODE'/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
C              CHECK IF OPTION EXISTS
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
      WRITE(CHMAIL,2000)IAT
      CALL GMAIL(0,0)
      GO TO 99
C              CHECK IF VOLUME EXISTS
C
  10  IF(JVOLUM.GT.0)GO TO 20
      WRITE(CHMAIL,1000)NAME
      CALL GMAIL(0,0)
      GO TO 99
C
  20  IF(NAME(1:1).EQ.'*')GO TO 40
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVOLU)
      IF(IVOLU.GT.0)GO TO 30
      WRITE(CHMAIL,1000)NAME
      CALL GMAIL(0,0)
      GO TO 99
C
  30  JVO=LQ(JVOLUM-IVOLU)
      NPAR=Q(JVO+IPN)
      ISVALS=Q(JVO+NPAR+IPS-1+IOP)
      GO TO 99
C
C             SET OPTION VALUE FOR ALL DEFINED VOLUMES
C
  40  CONTINUE
*      DO 50 I=1,NVOLUM
      I=2
      JVO=LQ(JVOLUM-I)
      NPAR=Q(JVO+IPN)
      ISVALS=Q(JVO+NPAR+IPS-1+IOP)
*  50  CONTINUE
C
 1000 FORMAT(' GSATT VOLUME ',A4,' DOES NOT EXIST ')
 2000 FORMAT(' GSATT OPTION ',A4,' NOT IMPLEMENTED ')
  99  RETURN
      END

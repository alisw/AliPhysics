*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDOPT(IOPTC,IVALC)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Set options for the drawing package                      *
C.    *                                                                *
C.    *       IOPTN = option (input)                                   *
C.    *       IVALN = value (input)                                    *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini  J.Salt  S.Giani                     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcspee.inc"
#include "geant321/gcrayt.inc"
#include "geant321/gcvdma.inc"
#include "geant321/gcpixe.inc"
C
      CHARACTER*4 IOPTC,IVALC
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH(IVALC,IVAL,4,4)
      IF (IOPTC.EQ.'THRZ') THEN
         IF (IVALC.EQ.'ON  ') THEN
            ITHRZ=IVAL
         ELSE IF (IVALC.EQ.'OFF ') THEN
            ITHRZ=IVAL
         ELSE IF (IVALC.EQ.'180 ') THEN
            ITHRZ=IVAL
         ELSE IF (IVALC.EQ.'360 ') THEN
            ITHRZ=IVAL
         ELSE
            WRITE (CHMAIL,10200) IOPTC,ITHRZ
            CALL GMAIL(0,0)
         ENDIF
      ELSE IF (IOPTC.EQ.'PROJ') THEN
         IF (IVALC.EQ.'PARA') THEN
            IPRJ=IVAL
         ELSE IF (IVALC.EQ.'PERS') THEN
            IPRJ=IVAL
         ELSE
            WRITE (CHMAIL,10200) IOPTC,IPRJ
            CALL GMAIL(0,0)
         ENDIF
      ELSE IF (IOPTC.EQ.'TRAK') THEN
         IF (IVALC.EQ.'LINE') THEN
            ITRKOP=IVAL
         ELSE IF (IVALC.EQ.'POIN') THEN
            ITRKOP=IVAL
         ELSE
            WRITE (CHMAIL,10200) IOPTC,ITRKOP
            CALL GMAIL(0,0)
         ENDIF
*JS
      ELSE IF (IOPTC.EQ.'HIDE') THEN
         IF (IVALC.EQ.'ON  ') THEN
            IHIDEN=IVAL
         ELSE IF (IVALC.EQ.'OFF ') THEN
            IHIDEN=IVAL
         ELSE
            WRITE (CHMAIL,10200) IOPTC,IHIDEN
            CALL GMAIL(0,0)
         ENDIF
*JS
**SG
      ELSE IF (IOPTC.EQ.'SHAD') THEN
         IF (IVALC.EQ.'ON  ') THEN
            LEP=10
         ELSE IF (IVALC.EQ.'OFF ') THEN
            LEP=1
         ENDIF
 
      ELSE IF (IOPTC.EQ.'EDGE') THEN
         IF (IVALC.EQ.'OFF ') THEN
            LLEP=ABS(LEP)
            IF(LLEP.EQ.10)LEP=11
         ELSE IF (IVALC.EQ.'ON  ') THEN
            LLEP=ABS(LEP)
            IF(LLEP.EQ.11)LEP=10
            IF(LLEP.EQ.1)LEP=1
         ENDIF
 
      ELSE IF (IOPTC.EQ.'RAYT') THEN
         IF (IVALC.EQ.'ON  ') THEN
             RAYTRA=1.
             IOMBRA=0
#if defined(CERNLIB_OLD)
             RAYTRA=0.
#endif
         ELSE IF (IVALC.EQ.'OFF ') THEN
             RAYTRA=0.
         ENDIF
 
      ELSE IF (IOPTC.EQ.'MAPP') THEN
         IF (IVALC.EQ.'0   ') THEN
             NMAP=0
         ELSE IF (IVALC.EQ.'1   ') THEN
             NMAP=1
         ELSE IF (IVALC.EQ.'2   ') THEN
             NMAP=2
         ELSE IF (IVALC.EQ.'3   ') THEN
             NMAP=3
         ELSE IF (IVALC.EQ.'4   ') THEN
             NMAP=4
         ENDIF
 
      ELSE IF (IOPTC.EQ.'USER') THEN
         IF (IVALC.EQ.'ON  ') THEN
             IMYSE=1
             CALL GSATT('*   ','SEEN',-10)
             CALL GSATT('*   ','COLO',-10)
             CALL GSATT('*   ','LSTY',-10)
         ELSE IF (IVALC.EQ.'OFF ') THEN
             IMYSE=0
             CALL GSATT('*   ','SEEN',1)
             CALL GSATT('*   ','COLO',1)
             CALL GSATT('*   ','LSTY',1)
         ENDIF
**SG
      ELSE
         WRITE (CHMAIL,10000) IOPTC
         CALL GMAIL(0,0)
         WRITE (CHMAIL,10100)
         CALL GMAIL(0,0)
      ENDIF
C
10000 FORMAT(' GDOPT: ',A4,' is not a valid option')
*JS
10100 FORMAT('    Valid options : THRZ , PROJ , TRAK , HIDE , SHAD')
*JS
10200 FORMAT(' GDOPT: Option ',A4,' is ',A4)
*SG
10300 FORMAT(' GDOPT: Option ',A4,' is ',A4)
10400 FORMAT(' Please, compute and set the right size of Zebra',
     +       ' store, needed for your drawing, with LHC option',
     +       ' OFF. If it is done, good luck !!!')
*SG
      END
 

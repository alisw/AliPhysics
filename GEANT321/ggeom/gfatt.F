*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFATT(NAME,IAT,IVAL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        Fetch value IVAL of attribute IAT                       *
C.    *        for the volume NAME                                     *
C.    *                                                                *
C.    *    ==>Called by : GDPRTR,GDSPEC,GDTREE,GDTR1                   *
C.    *       Author  R.Brun, P.Zanarini  **********                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 NAME,IAT,KATT(10)
      DIMENSION LATT(10)
      SAVE IFIRST,KATT,LATT
      DATA KATT/'WORK','SEEN','LSTY','LWID','COLO','FILL','SET ',
     +'DET ','DTYP','NODE'/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
C              CHECK IF OPTION EXISTS
C
      IF(NAME.EQ.' ')GO TO 99
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH(KATT,LATT,4,40)
      ENDIF
C
      CALL GLOOK(IAT,LATT,10,IOP)
      IF(IOP.GT.0)GO TO 10
      WRITE(CHMAIL,2000)IAT
      CALL GMAIL(0,0)
      GO TO 99
C
C              CHECK IF VOLUME EXISTS
C
  10  CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVOLU)
      IF(IVOLU.GT.0)GO TO 30
      WRITE(CHMAIL,1000)NAME
      CALL GMAIL(0,0)
      GO TO 99
C
  30  JVO=LQ(JVOLUM-IVOLU)
      NPAR=Q(JVO+5)
      IVAL=Q(JVO+NPAR+6+IOP)
      GO TO 99
C
 1000 FORMAT(' GFATT VOLUME ',A4,' DOES NOT EXIST ')
 2000 FORMAT(' GFATT OPTION ',A4,' NOT IMPLEMENTED ')
  99  RETURN
      END

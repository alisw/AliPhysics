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
      SUBROUTINE GSCHK(NAME,KSHAPE,NPAR,ISHAPE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *           CHECK SHAPE PARAMETERS                               *
C.    *           RETURN SHAPE NUMBER CORRESPONDING TO 'KSHAPE'        *
C.    *                                                                *
C.    *    ==>Called by : GSVOLU                                       *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
      CHARACTER*4 NAME,KSHAP(30)
      CHARACTER*4    KSHAPE
      DIMENSION NPARS(30)
      SAVE KSHAP,NPARS
      DATA KSHAP/'BOX ','TRD1','TRD2','TRAP','TUBE','TUBS','CONE',
     +           'CONS','SPHE','PARA','PGON','PCON','ELTU','HYPE',
     + 13*'    ','GTRA','CTUB','    '/
      DATA NPARS/3,4,5,11,3,5,5,7,6,6,10,9,3,4,13*0,12,11,0/
C.
C.    ------------------------------------------------------------------
C.
      ISHAPE = 0
      DO 10 JSH=1, 30
         IF(KSHAPE.EQ.KSHAP(JSH)) THEN
            ISHAPE=JSH
            GOTO 20
         ENDIF
   10 CONTINUE
C
      WRITE(CHMAIL,10000)KSHAPE,NAME
      CALL GMAIL(1,0)
      GO TO 999
C
C              CHECK PARAMETERS
C
   20 IF(NPAR.LT.NPARS(ISHAPE).AND.NPAR.NE.0) THEN
C
         WRITE(CHMAIL,10100)KSHAPE,NAME,NPAR
         CALL GMAIL(1,0)
         ISHAPE=-1
         GO TO 999
      ENDIF
C
10000 FORMAT(' ***** SHAPE TYPE ',A4,' VOLUME ',A4,
     +       ' NOT VALID *****')
10100 FORMAT(' ***** SHAPE ',A4,' VOLUME ',A4,
     +       ' INCONSISTENT NUMBER OF PARAMETERS ',I5,'*****')
  999 END

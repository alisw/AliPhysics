*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 13/12/94  15.26.17  by  S.Giani
*-- Author :
      SUBROUTINE GPVOLU(NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print VOLUMES data structures JVOLUM          *
C.    *       NUMB     Volume number                                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      COMMON/FMOTH/INGLOB,IVOMGL
      CHARACTER*4 ISHAP(30),NAME,MOTHER
      DIMENSION PAR(100),ATT(20)
      SAVE ISHAP
      DATA ISHAP/'BOX ','TRD1','TRD2','TRAP','TUBE','TUBS','CONE',
     +           'CONS','SPHE','PARA','PGON','PCON','ELTU','HYPE',
     +13*'    ', 'GTRA','CTUB','    '/
C.
C.    ------------------------------------------------------------------
      IF (JVOLUM.LE.0) GO TO 99
      IF (NUMB  .GT.0) GO TO 10
C
      WRITE (CHMAIL,1000)
      CALL GMAIL(0,0)
      WRITE (CHMAIL,1001)
      CALL GMAIL(0,1)
C
      N1     = 1
      N2     = NVOLUM
      GO TO 20
C
   10 WRITE (CHMAIL,1001)
      CALL GMAIL(0,1)
C
      N1     = NUMB
      N2     = NUMB
C
   20 DO 30 I=N1,N2
      JVO = LQ(JVOLUM-I)
      NIN=Q(JVO+3)
      IF (JVO.LE.0)   GO TO 30
C
      CALL UHTOC(IQ(JVOLUM+I),4,NAME,4)
      CALL GFMOTH(NAME,MOTHER,KONLY)
      JVOMOT=LQ(JVOLUM-IVOMGL)
      MNIN=Q(JVOMOT+3)
      IF(MNIN.LE.0)THEN
        NMBR=1
      ELSE
        JIN = LQ(JVOMOT-INGLOB)
        NMBR=Q(JIN+3)
      ENDIF
 
      IS   = Q(JVO+2)
      NMED = Q(JVO+4)
C
      CALL GFPARA(NAME,NMBR,1,NPAR,NATT,PAR,ATT)
      IEND10=10
      IF(NPAR.LT.IEND10)IEND10=NPAR
      WRITE(CHMAIL,1002)I,NAME,NMED,ISHAP(IS),NPAR,(PAR(J),J=1,IEND10)
      CALL GMAIL(0,0)
      DO 25 I10=11,NPAR,10
        IEND10=I10+9
        IF (NPAR.LT.IEND10) IEND10=NPAR
        WRITE (CHMAIL,1003) (PAR(J),J = I10,IEND10)
        CALL GMAIL(0,0)
   25 CONTINUE
   30 CONTINUE
C
 1000 FORMAT ('0',51('='),5X,' VOLUMES ',6X,50('='))
 1001 FORMAT ('0','VOLUME NAME  NUMED SHAPE NPAR  PARAMETERS')
 1002 FORMAT (' ',I6,1X,A4,2X,I3,3X,A4,I5,2X,10E10.3)
 1003 FORMAT (31X,10E10.3)
  99  END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GEDITV(IMENU)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Edit volumes (only for interactive version)              *
C.    *                                                                *
C.    *       IMENU = option selected from menu (input)                *
C.    *                                                                *
C.    *    ==>Called by : GINC3                                        *
C.    *       Author    P.Zanarini  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      CHARACTER*4 CHNAME,CHNEW,NAMDIV,NAMMOT
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH('    ',IBLA,4,4)
C
      CALL KUPROC('Give volume NAME',CHNAME,NCH)
      CALL UCTOH(CHNAME,NAME,4,NCH)
      IVO=IUCOMP(NAME,IQ(JVOLUM+1),NVOLUM)
      IF (IVO.LE.0) GO TO 999
      JVO=LQ(JVOLUM-IVO)
C
C             Get IVOMOT,NIN,JIN
C
      IF (IMENU.GE.4.AND.IMENU.LE.6) THEN
         CALL KUPROI('Give copy NR',NR)
         DO 20 IVOMOT=1,NVOLUM
            JVOMOT=LQ(JVOLUM-IVOMOT)
            NIN=Q(JVOMOT+3)
            IF (NIN.LT.0) THEN
               NIN=1
               IDIV=1
            ELSE
               IDIV=0
            ENDIF
            DO 10 IN=1,NIN
               JIN=LQ(JVOMOT-IN)
               IVOSON=Q(JIN+2)
               IF (IDIV.EQ.1) THEN
                  IF (IVOSON.EQ.IVO) GO TO 30
               ELSE
                  NRSON=Q(JIN+3)
                  IF (IVOSON.EQ.IVO.AND.NRSON.EQ.NR) GO TO 30
               ENDIF
   10       CONTINUE
   20    CONTINUE
         GO TO 999
   30    CONTINUE
      ENDIF
C
C             Get IVOMOT,JDIV
C
      IF (IMENU.GE.7.AND.IMENU.LE.8) THEN
         DO 40 IVOMOT=1,NVOLUM
            JVOMOT=LQ(JVOLUM-IVOMOT)
            NIN=Q(JVOMOT+3)
            IF (NIN.GE.0) GO TO 40
            JDIV=LQ(JVOMOT-1)
            IVOSON=Q(JDIV+2)
            IF (IVOSON.EQ.IVO) GO TO 50
   40    CONTINUE
         GO TO 999
   50    CONTINUE
      ENDIF
C
      IF (IMENU.EQ.1) THEN
C
C             Modify shape parameters PAR given by GSVOLU
C
         NP=Q(JVO+5)
         DO 60 I=1,NP
            PAR=Q(JVO+6+I)
            WRITE (CHMAIL,1100) I,PAR
            CALL GMAIL(0,0)
            CALL KUPROR('Give new value',PAR)
            Q(JVO+6+I)=PAR
   60    CONTINUE
 
C
      ELSE IF (IMENU.EQ.2) THEN
C
C             Modify NAME given by GSVOLU
C
         CALL KUPROC('Give new NAME',CHNEW,NCH)
         NEWNAM=IBLA
         CALL UCTOH(CHNEW,NEWNAM,4,NCH)
         IQ(JVOLUM+IVO)=NEWNAM
C
      ELSE IF (IMENU.EQ.3) THEN
C
C             Delete NAME given by GSVOLU
C
   70    CONTINUE
         DO 90 IVOMOT=1,NVOLUM
            JVOMOT=LQ(JVOLUM-IVOMOT)
            NIN=Q(JVOMOT+3)
            IF (NIN.LT.0) NIN=1
            DO 80 IN=1,NIN
               JIN=LQ(JVOMOT-IN)
               IVOSON=Q(JIN+2)
               NR=Q(JIN+3)
               IF (IVOSON.EQ.IVO) GO TO 100
   80       CONTINUE
   90    CONTINUE
         GO TO 110
C
  100    CONTINUE
C
C             Unlink NAME,NR
C
#if !defined(CERNLIB_USRJMP)
         CALL GUNLIV(IVO,NR,IVOMOT)
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT3(JUNLIV,IVO,NR,IVOMOT)
#endif
         WRITE (CHMAIL,1000) NAME
         CALL GMAIL(0,0)
 1000   FORMAT (' *** GEDITV: ',A4,' UNLINKED')
C
C             Try another link
C
         GO TO 70
C
  110    CONTINUE
C
C             No more links; now delete NAME
C
         IQ(JVOLUM+IVO)=IBLA
C
      ELSE IF (IMENU.EQ.4) THEN
C
C             Unlink NAME,NR given by GSPOS/GSDIV
C
#if !defined(CERNLIB_USRJMP)
         CALL GUNLIV(IVO,NR,IVOMOT)
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT3(JUNLIV,IVO,NR,IVOMOT)
#endif
C
      ELSE IF (IMENU.EQ.5) THEN
C
C             Modify X0,Y0,Z0 of NAME,NR given by GSPOS
C
         X0=Q(JIN+5)
         Y0=Q(JIN+6)
         Z0=Q(JIN+7)
         CALL KUPROR('Give X0',X0)
         CALL KUPROR('Give Y0',Y0)
         CALL KUPROR('Give Z0',Z0)
         Q(JIN+5)=X0
         Q(JIN+6)=Y0
         Q(JIN+7)=Z0
C
      ELSE IF (IMENU.EQ.6) THEN
C
C             Modify IROT of NAME,NR given by GSPOS
C
         IROT=Q(JIN+4)
         CALL KUPROI('Give IROT',IROT)
         Q(JIN+4)=IROT
C
      ELSE IF (IMENU.EQ.7.OR.IMENU.EQ.8) THEN
C
         IF (IMENU.EQ.7) THEN
C
C             Modify NDIV given by GSDIV
C
            NDIV=Q(JDIV+3)
            CALL KUPROI('Give NDIV',NDIV)
            Q(JDIV+3)=NDIV
C
         ELSE
C
C             Modify IAXIS given by GSDIV
C
            IAXIS=Q(JDIV+1)
            CALL KUPROI('Give IAXIS',IAXIS)
            Q(JDIV+1)=IAXIS
C
         ENDIF
C
C             Unlink and delete NAME
C
         Q(JVOMOT+3)=0
         CALL MZDROP(IXCONS,LQ(JVOMOT-1),' ')
         JV = LQ(JVOLUM-IVOMOT)
         CALL MZPUSH(IXCONS,JV,-1,0,'I')
         CALL UHTOC(IQ(JVOLUM+IVO),4,NAMDIV,4)
         IQ(JVOLUM+IVO)=IBLA
C
C             Redivide (division is now at NVOLUM-th position)
C
         CALL UHTOC(IQ(JVOLUM+IVOMOT),4,NAMMOT,4)
         NDIV=Q(JDIV+3)
         IAXIS=Q(JDIV+1)
         CALL GSDVN(NAMDIV,NAMMOT,NDIV,IAXIS)
C
C             Swap new division with old one (links + names)
C
         CALL DZSWAP(IXCONS,LQ(JVOLUM-NVOLUM),LQ(JVOLUM-IVO),' ')
         IQ(JVOLUM+IVO)=IQ(JVOLUM+NVOLUM)
         IQ(JVOLUM+NVOLUM)=IBLA
         JVOMOT=LQ(JVOLUM-IVOMOT)
         JDIV=LQ(JVOMOT-1)
         Q(JDIV+2)=IVO
         CALL UCTOH(NAMDIV,IQ(JVOLUM+IVO),4,4)
C
C             Delete definitely old division
C
         CALL MZDROP(IXCONS,LQ(JVOLUM-NVOLUM),' ')
         CALL MZPUSH(IXCONS,JVOLUM,-1,-1,'I')
         NVOLUM=NVOLUM-1
C
      ENDIF
C
 1100 FORMAT('  PAR(',I2,') =',F10.3)
  999 RETURN
      END

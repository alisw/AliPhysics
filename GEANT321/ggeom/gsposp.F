*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSPOSP(NAME,NR,MOTHER,X,Y,Z,IROT,KONLY,UPAR,NP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Place a copy of generic volume 'NAME' with user number    *
C.    *      'NR' inside 'MOTHER', with its parameters UPAR(1..NP)     *
C.    *                                                                *
C.    *          JVO=pointer to mother volume                          *
C.    *          JIN=pointer to the copy 'NAME','NR'                   *
C.    *          JIN=LQ(JVO-IN)                                        *
C.    *                                                                *
C.    *             Q(JIN+1)=NENTRY                                    *
C.    *             Q(JIN+2)=VOLUME NUMBER                             *
C.    *             Q(JIN+3)=USER NUMBER                               *
C.    *             Q(JIN+4)=IROT                                      *
C.    *             Q(JIN+5)=X                                         *
C.    *             Q(JIN+6)=Y                                         *
C.    *             Q(JIN+7)=Z                                         *
C.    *             Q(JIN+8)=ONLY                                      *
C.    *             Q(JIN+9)=NPAR                                      *
C.    *             Q(JIN+10 ..)=PAR ..                                *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *         Authors R.Brun, F.Bruyant,  A.McPherson  *********     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconsp.inc"
      CHARACTER*4 NAME,MOTHER,KONLY
      DIMENSION UPAR(*),PAR(100)
C.
C.    ------------------------------------------------------------------
C.
C              Check if volume master bank exists
C
      IF(JVOLUM.GT.0)GO TO 10
      WRITE(CHMAIL,1000)
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if mother volume exists
C
  10  CALL GLOOK(MOTHER,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.GT.0)GO TO 20
      WRITE(CHMAIL,2000)MOTHER
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if NAME volume exists
C
  20  CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IN)
      IF(IN.GT.0)GO TO 30
      WRITE(CHMAIL,2000)NAME
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if rotation matrix exists
C
  30  IF(IROT.LE.0)GO TO 50
      IF(JROTM.GT.0)GO TO 40
      WRITE(CHMAIL,3000)IROT
      CALL GMAIL(0,0)
      GO TO 99
  40  IF(LQ(JROTM-IROT).GT.0)GO TO 50
      WRITE(CHMAIL,3000)IROT
      CALL GMAIL(0,0)
      GO TO 99
C
C              Check if mother is not divided
C
  50  JIN=LQ(JVOLUM-IN)
      ISH=Q(JIN+2)
      JVO=LQ(JVOLUM-IVO)
      ICOPY=1
      NIN=Q(JVO+3)
      IF(NIN.GE.0)GO TO 60
      WRITE(CHMAIL,4000)MOTHER
      CALL GMAIL(0,0)
      GO TO 99
*
* *** Copy user parameters into local array PAR
  60  NPAR=NP
      IF (ISH.EQ. 4) NPAR=35
      IF (ISH.EQ.28) NPAR=30
      CALL UCOPY(UPAR,PAR,NP)
*
* *** Check if ('NAME',NUMBER') exists
      IF(NIN.EQ.0)GO TO 80
      DO 70 I=1,NIN
      JIN=LQ(JVO-I)
      IF(Q(JIN+2).NE.IN)GO TO 70
      IF(Q(JIN+3).NE.NR)GO TO 70
      WRITE(CHMAIL,5000)NAME,NR
      CALL GMAIL(0,0)
      GO TO 90
  70  CONTINUE
      ICOPY=NIN+1
*
* *** Create bank for that copy
  80  NINL=IQ(JVO-2)
      IF(ICOPY.GT.NINL)CALL MZPUSH(IXCONS,JVO,50,0,'I')
      CALL MZBOOK(IXCONS,JIN,JVO,-ICOPY,'VOPP',1,1,NPAR+9,3,0)
      IF(IEOTRI.NE.0)GO TO 95
      IQ(JIN-5)=100*IVO+ICOPY
      Q(JVO+3)=Q(JVO+3)+1
*
* *** Now store parameters into bank area
  90  Q(JIN+2)=IN
      Q(JIN+3)=NR
      Q(JIN+4)=IROT
      Q(JIN+5)=X
      Q(JIN+6)=Y
      Q(JIN+7)=Z
      IF(KONLY.EQ.'ONLY')Q(JIN+8)=1.
      Q(JIN+9) = NPAR
*
      IF (ISH.EQ.4) THEN
*        Trapezoid
         TTH= TAN(PAR(2)*DEGRAD)
         PHI    = PAR(3)*DEGRAD
         PAR(2) = TTH*COS(PHI)
         PAR(3) = TTH*SIN(PHI)
         PAR(7) = TAN(PAR(7) *DEGRAD)
         PAR(11)= TAN(PAR(11)*DEGRAD)
         CALL GNOTR1 (PAR)
      ELSE IF (ISH.EQ.10) THEN
*        Parallelepiped.
         PAR(4)=TAN(PAR(4)*DEGRAD)
         TTH=TAN(PAR(5)*DEGRAD)
         PH=PAR(6)*DEGRAD
         PAR(5)=TTH*COS(PH)
         PAR(6)=TTH*SIN(PH)
      ELSE IF (ISH.EQ.28) THEN
*        General twisted trapezoid.
         CALL GTRAIN(UPAR,PAR)
      ENDIF
*
      CALL UCOPY(PAR,Q(JIN+10),NPAR)
      GO TO 99
*
*     Not enough space
  95  WRITE(CHMAIL,6000)NAME,NR,MOTHER
      CALL GMAIL(0,0)
C
 1000 FORMAT(' ***** GSPOSP CALLED AND NO VOLUMES DEFINED *****')
 2000 FORMAT(' ***** GSPOSP VOLUME ',A4,' DOES NOT EXISTS *****')
 3000 FORMAT(' ***** GSPOSP ROTATION MATRIX',I5,' DOES NOT EXIST *****')
 4000 FORMAT(' ***** GSPOSP MOTHER ',A4,' ALREADY DIVIDED *****')
 5000 FORMAT(' ***** GSPOSP COPY ',A4,' NUMBER ',I5,
     +       ' ALREADY CREATED IN ',A4,' *****')
 6000 FORMAT(' ***** GSPOSP NOT ENOUGH SPACE TO STORE COPY ',A4,
     +       ' NUMBER ',I5,' IN ',A4,' *****')
  99  RETURN
      END

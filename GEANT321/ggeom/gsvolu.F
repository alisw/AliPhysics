*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSVOLU(KNAME,JSHAPE,NMED,UPAR,NP,IVOLU)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        CREATES A NEW VOLUME                                    *
C.    *                                                                *
C.    *          JVO=LQ(JVOLUM-IVOLU)                                  *
C.    *                                                                *
C.    *            Q(JVO+1)=ISEARC (SET TO 0 BY DEFAULT)               *
C.    *            Q(JVO+2)=ISHAPE                                     *
C.    *            Q(JVO+3)=NIN                                        *
C.    *            Q(JVO+4)=NMED                                       *
C.    *            Q(JVO+5)=NPAR                                       *
C.    *            Q(JVO+6)=NATT                                       *
C.    *            Q(JVO+7)=PAR.....                                   *
C.    *            Q(JVO+7+NPAR)=ATT.....                              *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 KNAME,JSHAPE
      DIMENSION UPAR(50),PAR(50),ATT(20)
      SAVE ATT,NATT
      DATA ATT/1.,1.,1.,1.,1.,15*0./
      DATA NATT/10/
C.
C.    ------------------------------------------------------------------
C.
C
C              Copy user parameters into local array PAR
C
      NPAR=NP
      IF (NP.GT.0) THEN
         IF(JSHAPE.EQ.'TRAP') NPAR=35
         IF(JSHAPE.EQ.'GTRA') NPAR=30
         CALL UCOPY(UPAR,PAR,NP)
      ENDIF
C
      IVOLU=0
C
C              CHECK SHAPE VALIDITY
C
      CALL GSCHK ( KNAME, JSHAPE, NPAR, ISHAPE )
C
      IF(ISHAPE.LE.0)GO TO 99
      IF(JVOLUM.GT.0)GO TO 10
C
C              CREATE THE MOTHER MEDIA BANK
C
      CALL MZBOOK(IXCONS,JVOLUM,JVOLUM,1,'VOLU',400,400,400,5,0)
      IVO=1
      NVOLUM=1
      GO TO 30
C
C              CHECK IF SUCH A VOLUME ALREADY DEFINED
C
  10  NVOL=IQ(JVOLUM-2)
      CALL GLOOK(KNAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.LE.0)GO TO 20
      IF(LQ(JVOLUM-IVO).GT.0) THEN
         WRITE(CHMAIL,10010) KNAME, IVO
10010 FORMAT(' **** GSVOLU: Redefinition of volume ',
     +       A4,' IVO = ',I6)
         CALL GMAIL(0,0)
         CALL MZDROP(IXCONS,LQ(JVOLUM-IVO),' ')
      ENDIF
      GO TO 30
C
  20  NVOLUM=NVOLUM+1
      IF(NVOLUM.GT.NVOL)CALL MZPUSH(IXCONS,JVOLUM,50,50,'I')
      IVO=NVOLUM
C
C              NOW CREATE THE VOLUME BANK
C
  30  CALL MZBOOK(IXCONS,JVO,JVOLUM,-IVO,'VOL1',50,50,9+NPAR+NATT,3,0)
      CALL UCTOH(KNAME,IQ(JVOLUM+IVO),4,4)
C
C              COPY PARAMETERS IN DATA AREA
C
      IVOLU=IVO
      Q(JVO+2)=ISHAPE
      Q(JVO+4)=NMED
      IF(NPAR.LE.0)GO TO 99
*
      IF (ISHAPE.EQ.4) THEN
*        Trapezoid
         TTH= TAN(PAR(2)*DEGRAD)
         PHI    = PAR(3)*DEGRAD
         PAR(2) = TTH*COS(PHI)
         PAR(3) = TTH*SIN(PHI)
         PAR(7) = TAN(PAR(7) *DEGRAD)
         PAR(11)= TAN(PAR(11)*DEGRAD)
         CALL GNOTR1 (PAR)
      ELSE IF (ISHAPE.EQ.10) THEN
*        Parallelepiped change angles to tangents.
         PAR(4)=TAN(PAR(4)*DEGRAD)
         TTH=TAN(PAR(5)*DEGRAD)
         PH=PAR(6)*DEGRAD
         PAR(5)=TTH*COS(PH)
         PAR(6)=TTH*SIN(PH)
      ELSE IF (ISHAPE.EQ.28) THEN
*        General twisted trapezoid.
         CALL GTRAIN(UPAR,PAR)
      ENDIF
*
  99  CONTINUE
*
      CALL GSIPAR(JVO,0,NPAR,NATT,PAR,ATT)
*
      END

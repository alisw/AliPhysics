*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GREUCL (LUN,FILNAM)
*
*
*     ******************************************************************
*     *                                                                *
*     *  Read in the geometry of the detector in EUCLID file format    *
*     *                                                                *
*     *                                                                *
*     *       Author : M. Maire                                        *
*     *                                                                *
*     ******************************************************************
*
      CHARACTER*(*) FILNAM
      CHARACTER*80  FILEXT
      CHARACTER    KEY*4, CARD*76
      CHARACTER*20 NATMED, NAMATE
      CHARACTER*4  NAME, MOTHER, SHAPE, KONLY
      REAL PAR(100)
*
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
*
      IT2=LNBLNK(FILNAM)
#if !defined(CERNLIB_IBM)
      FILEXT=FILNAM
#endif
#if defined(CERNLIB_IBM)
      IT1=INDEX(FILNAM,'.')
      IF (IT1.GT.0) THEN
         FILEXT='/'//FILNAM(1:IT1-1)//' '//FILNAM(IT1+1:IT2)//' A1'
         CALL CLTOU(FILEXT)
      ELSE
         WRITE (CHMAIL,10000) FILNAM
         CALL GMAIL (0,0)
         RETURN
      ENDIF
#endif
*
      OPEN (UNIT=LUN,FILE=FILEXT,STATUS='UNKNOWN',FORM='FORMATTED')
*
   10 READ (LUN,'(A4,A76)',END=20,ERR=20) KEY,CARD
*
      IF (KEY.EQ.'TMED') THEN
         READ (CARD,*)    ITMED, NATMED
         CALL GCKMAT (ITMED, NATMED)
*
      ELSE IF (KEY.EQ.'ROTM') THEN
         READ (CARD,*) IROT, TETA1, PHI1, TETA2, PHI2, TETA3, PHI3
         CALL GSROTM ( IROT, TETA1, PHI1, TETA2, PHI2, TETA3, PHI3 )
*
      ELSE IF (KEY.EQ.'VOLU') THEN
         READ (CARD,*) NAME, SHAPE, NUMED, NPAR
         IF (NPAR.GT.0) READ (LUN, *) (PAR(I),I=1,NPAR)
         CALL GSVOLU ( NAME, SHAPE, NUMED, PAR, NPAR, IVOLU)
*
      ELSE IF (KEY.EQ.'DIVN') THEN
         READ (CARD,*) NAME, MOTHER, NDIV, IAXE
         CALL GSDVN  ( NAME, MOTHER, NDIV, IAXE )
*
      ELSE IF (KEY.EQ.'DVN2') THEN
         READ (CARD,*) NAME, MOTHER, NDIV, IAXE, ORIG, NUMED
         CALL GSDVN2 ( NAME, MOTHER, NDIV, IAXE, ORIG, NUMED )
*
      ELSE IF (KEY.EQ.'DIVT') THEN
         READ (CARD,*) NAME, MOTHER, STEP, IAXE, NUMED, NDVMX
         CALL GSDVT  ( NAME, MOTHER, STEP, IAXE, NUMED, NDVMX )
*
      ELSE IF (KEY.EQ.'DVT2') THEN
         READ (CARD,*) NAME, MOTHER, STEP, IAXE, ORIG, NUMED, NDVMX
         CALL GSDVT2 ( NAME, MOTHER, STEP, IAXE, ORIG, NUMED, NDVMX )
*
      ELSE IF (KEY.EQ.'POSI') THEN
         READ (CARD,*) NAME, NR, MOTHER, XO, YO, ZO, IROT, KONLY
         CALL GSPOS  ( NAME, NR, MOTHER, XO, YO, ZO, IROT, KONLY )
*
      ELSE IF (KEY.EQ.'POSP') THEN
         READ (CARD,*) NAME, NR, MOTHER, XO, YO, ZO, IROT, KONLY, NPAR
         IF (NPAR.GT.0) READ (LUN, *) (PAR(I),I=1,NPAR)
         CALL GSPOSP ( NAME, NR, MOTHER, XO,YO,ZO,IROT,KONLY,PAR,NPAR)
      ENDIF
*
      IF (KEY.NE.'END') GO TO 10
      CLOSE (LUN)
*
      WRITE (CHMAIL,10100) FILNAM(1:IT2)
      CALL GMAIL (1,1)
*
      GOTO 999
*
   20 CONTINUE
      WRITE (CHMAIL,10200)
      CALL GMAIL (1,1)
      CLOSE (LUN)
*
10000 FORMAT (' *** GREUCL *** ',A,' is not a valib IBM name')
10100 FORMAT (' *** GREUCL *** file: ',A,' is now read in')
10200 FORMAT (' *** GREUCL *** reading error or premature end of file')
*
  999 END

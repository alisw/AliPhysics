*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXLIST
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Geometry commands                                         *
C.    *                                                                *
C.    *       Authors:   F.Carminati **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gclist.inc"
#include "geant321/gcrz.inc"
      PARAMETER (NUMLIS=11)
      DIMENSION LLIST1(20,9), NLIST1(9)
      DIMENSION LLIST2(20,2), NLIST2(2)
      EQUIVALENCE (NLIST1,NHSTA),(LLIST1,LHSTA)
      EQUIVALENCE (NLIST2,NRGET),(LLIST2,LRGET)
      CHARACTER*4  CHPAR, CHPATL, CHLNAM(NUMLIS)
      DATA CHLNAM / 'HSTA','GET ','SAVE','SETS','PRIN','GEOM',
     +              'VIEW','PLOT','STAT','RGET','RSAV'/
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      ILIST = 0
      DO 10 KLIST=1,NUMLIS
         IF(CHPATL.EQ.CHLNAM(KLIST)) THEN
            ILIST = KLIST
         ENDIF
   10 CONTINUE
      IF(ILIST.GT.9) THEN
         ILIST=ILIST-9
         IF(NPAR.GT.0) THEN
            CALL KUGETC(CHPAR,NCH)
            IF(CHPAR.EQ.'.') THEN
               CALL VBLANK(LLIST2(1,ILIST),20)
               NLIST2(ILIST) = 0
            ELSE
               NLIST2(ILIST) = NPAR
               CALL UCTOH(CHPAR,LLIST2(1,ILIST),4,4)
               DO 20 IPAR=2,NPAR
                  CALL KUGETC(CHPAR,NCH)
                  CALL UCTOH(CHPAR,LLIST2(IPAR,ILIST),4,4)
   20          CONTINUE
            ENDIF
         ELSE
            NDONE = MIN(NLIST2(ILIST),10)
            WRITE(CHMAIL,10000) CHPATL, NDONE, (LLIST2(IPAR,ILIST),
     +      IPAR= 1,NDONE)
            CALL GMAIL(1,0)
            IF(NLIST2(ILIST).GT.10) THEN
               WRITE(CHMAIL,10100) CHPATL, NLIST2(ILIST), (LLIST2(IPAR,
     +         ILIST), IPAR=11,NLIST2(ILIST))
               CALL GMAIL(0,0)
            ENDIF
         ENDIF
      ELSEIF(ILIST.GT.0) THEN
         IF(NPAR.GT.0) THEN
            CALL KUGETC(CHPAR,NCH)
            IF(CHPAR.EQ.'.') THEN
               CALL VBLANK(LLIST1(1,ILIST),20)
               NLIST1(ILIST) = 0
            ELSE
               NLIST1(ILIST) = NPAR
               CALL UCTOH(CHPAR,LLIST1(1,ILIST),4,4)
               DO 30 IPAR=2,NPAR
                  CALL KUGETC(CHPAR,NCH)
                  CALL UCTOH(CHPAR,LLIST1(IPAR,ILIST),4,4)
   30          CONTINUE
            ENDIF
         ELSE
            NDONE = MIN(NLIST1(ILIST),10)
            WRITE(CHMAIL,10000) CHPATL, NDONE, (LLIST1(IPAR,ILIST),
     +      IPAR= 1,NDONE)
            CALL GMAIL(1,0)
            IF(NLIST1(ILIST).GT.10) THEN
               WRITE(CHMAIL,10100) CHPATL, NLIST1(ILIST), (LLIST1(IPAR,
     +         ILIST), IPAR=11,NLIST1(ILIST))
               CALL GMAIL(0,0)
            ENDIF
         ENDIF
      ENDIF
*
10000 FORMAT(1X,A4,' list  1-',I2,10(2X,A4))
10100 FORMAT(1X,A4,' list 11-',I2,10(2X,A4))
      END

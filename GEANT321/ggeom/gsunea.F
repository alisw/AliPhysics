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
      SUBROUTINE GSUNEA(CHNAME,ISEARC)
C.
C.    *****************************************************************
C.    *                                                               *
C.    *    ROUTINE TO SPECIFY USER ORDERING FOR THE CONTENTS          *
C.    *    OF A GIVEN VOLUME, BY SETTING THE SEARCH FLAG >0           *
C.    *                                                               *
C.    *    ==>Called by : USER                                        *
C.    *         Author  R.Brun, F.Bruyant  *********                  *
C.    *                                                               *
C.    *****************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      CHARACTER*4 CHNAME
C.
C.             -----------------------------------------
C.
C            Check JVOLUM
*
      IF(JVOLUM.LE.0) THEN
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
*
      ELSE
*
*           find volume CHNAME.
         CALL GLOOK(CHNAME,IQ(JVOLUM+1),NVOLUM,IVO)
         IF(IVO.LE.0) THEN
            WRITE(CHMAIL,10100) CHNAME
            CALL GMAIL(0,0)
*
         ELSE
*
*          Set new value of ISEARC.
            IF (ISEARC.LE.0) THEN
               WRITE (CHMAIL, 10200) ISEARC
               CALL GMAIL (0, 0)
            ELSE
               JVO=LQ(JVOLUM-IVO)
               Q(JVO+1)=ISEARC
            ENDIF
         ENDIF
      ENDIF
*
10000 FORMAT(' ***** GSUNEA - No volumes defined')
10100 FORMAT(' ***** GSUNEA - Volume ',A4,' does not exist')
10200 FORMAT(' ***** GSUNEA - ISEARC = ',I10,' should be positive')
      END

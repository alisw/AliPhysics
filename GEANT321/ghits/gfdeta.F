*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFDETA (IUSET,IUALI,IALI)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Returns the position IALI of the detector alias IUALI      *
C.    *                                                                *
C.    *    ==>Called by : <USER>,                                      *
C.    *       Author    F.Bruyant  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 IUSET,IUALI
C
C.    ------------------------------------------------------------------
C.
      IF (JSET.LE.0) GO TO 90
      NSET = IQ(JSET-1)
      IF (NSET.LE.0) GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF (ISET.EQ.0) GO TO 90
      JS = LQ(JSET-ISET)
      NDET = IQ(JS-1)
      IF (NDET.LE.0) GO TO 90
      CALL GLOOK(IUALI,IQ(JS+1),NDET,IALI)
      IF (IALI.EQ.0) GO TO 95
      GO TO 99
C
   90 WRITE (CHMAIL, 1000) IUSET
      CALL GMAIL(0,0)
 1000 FORMAT (' ***** GFDETA ERROR FOR SET ',A4)
      GO TO 99
   95 WRITE (CHMAIL, 2000) IUSET,IUALI
      CALL GMAIL(0,0)
 2000 FORMAT (' ***** GFDETA ERROR FOR SET ',A4,' ALIAS ',A4,
     + ' NOT FOUND')
C
   99 RETURN
      END

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
      SUBROUTINE GFDETU(IUSET,IUDET,NUPAR,NW,UPAR)
C.
C.    *******************************************************************
C.    *                                                                 *
C.    *  Return in UPAR the first NUPAR user parameters of detector     *
C.    *  IUDET. NW is the total number of parameters(output), or zero   *
C.    *  if an  error occured.                                          *
C.    *                                                                 *
C.    *    ==>Called by : <USER>                                        *
C.    *       Author    W.Gebel  *********                              *
C.    *                                                                 *
C.    *******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION UPAR(1)
      CHARACTER*4 IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
      NW=0
      NSET=IQ(JSET-1)
      IF(NSET.LE.0) GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,IS)
      IF(IS.LE.0) GO TO 90
      JS=LQ(JSET-IS)
      ND=IQ(JS-1)
      IF(ND.LE.0) GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),ND,ID)
      IF(ID.LE.0) GO TO 95
      JD=LQ(JS-ID)
      JU=LQ(JD-3)
      IF(JU.NE.0)THEN
         NW=IQ(JU-1)
         IF(NUPAR.GT.0) CALL UCOPY(Q(JU+1),UPAR(1),NUPAR)
      ENDIF
      GO TO 99
C
C              Error
C
   90 WRITE (CHMAIL, 1000) IUSET
      CALL GMAIL(0,0)
 1000 FORMAT (' ***** GFDETU ERROR SET ',A4,' NOT FOUND')
      GO TO 99
   95 WRITE (CHMAIL, 2000) IUSET,IUDET
      CALL GMAIL(0,0)
 2000 FORMAT (' ***** GFDETU ERROR FOR SET ',A4,
     + ' DETECTOR ',A4,' NOT FOUND')
C
  99  RETURN
      END

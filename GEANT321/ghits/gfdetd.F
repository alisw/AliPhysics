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
      SUBROUTINE GFDETD(IUSET,IUDET,ND,NAMESD,NBITSD)
C.
C.    ******************************************************************
C.    *                                                                *
C.    * Returns digitisation  parameters for  detector IUDET  of set   *
C.    * IUSET.                                                         *
C.    *                                                                *
C.    *             Input parameters                                   *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user detector identifier                             *
C.    *                                                                *
C.    *             Output parameters                                  *
C.    * ND        number of elements per digitisation                  *
C.    * NAMESD    the  ND   variable  names  for   the  digitisation   *
C.    *           elements                                             *
C.    * NBITSD    the  ND  bit  numbers  for  packing  the  variable   *
C.    *           values.                                              *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION NBITSD(1)
      CHARACTER*4 IUSET,IUDET,NAMESD(1)
C.
C.    ------------------------------------------------------------------
C.
      IF(JSET.LE.0)GO TO 90
      NSET=IQ(JSET-1)
      IF(NSET.LE.0)GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF(ISET.LE.0)GO TO 90
      JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.LE.0)GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF(IDET.LE.0)GO TO 95
      JD=LQ(JS-IDET)
      JDD=LQ(JD-2)
      ND=IQ(JD+6)
C
      IF(ND.GT.0)THEN
         DO 10 I=1,ND
            CALL UHTOC(IQ(JDD+2*I-1),4,NAMESD(I),4)
            NBITSD(I)=IQ(JDD+2*I )
  10     CONTINUE
      ENDIF
      GO TO 99
C
C              Error
C
   90 WRITE (CHMAIL, 1000) IUSET
      CALL GMAIL(0,0)
 1000 FORMAT (' ***** GFDETD ERROR SET ',A4,' NOT FOUND')
      GO TO 99
   95 WRITE (CHMAIL, 2000) IUSET,IUDET
      CALL GMAIL(0,0)
 2000 FORMAT (' ***** GFDETD ERROR FOR SET ',A4,
     + ' DETECTOR ',A4,' NOT FOUND')
C
  99  RETURN
      END

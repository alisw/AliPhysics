*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GPSTAT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Print volume statistics                                  *
C.    *                                                                *
C.    *    ==>Called by : GLAST                                        *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
C.
C.      ---------------------------------------------------------------
C.
      IF(JGSTAT.LE.0)GO TO 999
      IF(NVOLUM.LE.0)GO TO 999
      WRITE(CHMAIL,10000)
      CALL GMAIL(3,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10200)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10100)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10000)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10300)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10400)
      CALL GMAIL(0,0)
      WRITE(CHMAIL,10000)
      CALL GMAIL(0,0)
C
      DO 10 I=1,NVOLUM
         JG=LQ(JGSTAT-I)
         IF(JG.GT.0) THEN
            WRITE(CHMAIL,10500)IQ(JVOLUM+I),(IQ(JG+J),J=1,9)
            CALL GMAIL(0,0)
         ENDIF
   10 CONTINUE
      WRITE(CHMAIL,10000)
      CALL GMAIL(0,2)
C
10000 FORMAT(1X,122('*'))
10100 FORMAT (' *',120X,'*')
10200 FORMAT (' *',43X,
     +'V O L U M E    S T A T I S T I C S',43X,'*')
10300 FORMAT (      ' *',3X,'VOLUME','   *',
     +3X,'      GINVOL      ',2X,'*',
     +3X,'      GMEDIA      ',2X,'*',
     +3X,'      GTNEXT      ',2X,'*',
     +3X,'      GTMEDI      ',2X,'*',
     +2X,'NSTEPS',3X,'*')
10400 FORMAT (      ' *',3X,' NAME ','   *',
     +3X,'TOTAL ','  *   ','NLEVEL',2X,'*',
     +3X,'TOTAL ','  *   ','NLEVEL',2X,'*',
     +3X,'TOTAL ','  *   ','NLEVEL',2X,'*',
     +3X,'TOTAL ','  *   ','NLEVEL',2X,'*',
     +2X,'NLEVEL',3X,'*')
10500 FORMAT(1X,'*',4X,A4,4X,9('*',I10,1X),'*')
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GPRUNG
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints Rung bank JRUNG                                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    M.Maire                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
*
      IF (JRUNG.NE.0) THEN
         NRUNG = IQ(JRUNG-1)
         WRITE (CHMAIL,1000) (IQ(JRUNG+I),I=1,10)
         CALL GMAIL(0,0)
         WRITE (CHMAIL,1001) (IQ(JRUNG+I),I=11,20)
         CALL GMAIL(0,0)
         WRITE (CHMAIL,1002) ( Q(JRUNG+I),I=21,30)
         CALL GMAIL(0,0)
         IF (LQ(JRUNG-1).NE.0) THEN
             WRITE(CHMAIL,2000)
             CALL GMAIL(0,0)
             JRU = LQ(JRUNG-1)
             NDT = IQ(JRU-1)
             NW  = MIN(NDT,9)
             WRITE(CHMAIL,3000) (Q(JRU+I),I=1,NW)
   15        CALL GMAIL(0,0)
             IF ((NDT-NW).GT.0) THEN
                NS = NW+1
                NW = MIN(NDT,NW+9)
                WRITE(CHMAIL,3000) (Q(JRU+I),I=NS,NW)
                GO TO 15
             ENDIF
         ENDIF
      ELSE
         WRITE(CHMAIL,4000)
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GPRUNG: IDRUN=',I4,2X,9I7)
 1001 FORMAT(20X,10I10)
 1002 FORMAT(20X,10F10.4)
 2000 FORMAT(14X,'User words')
 3000 FORMAT( 1X,9(2X,E12.5))
 4000 FORMAT(' **** GPRUNG: No Rung bank existing in memory ')
*
      END

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
      SUBROUTINE GPHEAD
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints Header bank JHEAD                                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    M.Maire                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
*
      IF (JHEAD.NE.0) THEN
         NHEAD = IQ(JHEAD-1)
         WRITE (CHMAIL,1000) (IQ(JHEAD+I),I=1,4)
         CALL GMAIL(0,0)
         WRITE (CHMAIL,1001) (IQ(JHEAD+I),I=5,NHEAD)
         CALL GMAIL(0,0)
         IF (LQ(JHEAD-1).NE.0) THEN
             WRITE(CHMAIL,2000)
             CALL GMAIL(0,0)
             JHU = LQ(JHEAD-1)
             NDT = IQ(JHU-1)
             NW  = MIN(NDT,9)
             WRITE(CHMAIL,3000) (Q(JHU+I),I=1,NW)
   15        CALL GMAIL(0,0)
             IF ((NDT-NW).GT.0) THEN
                NS = NW+1
                NW = MIN(NDT,NW+9)
                WRITE(CHMAIL,3000) (Q(JHU+I),I=NS,NW)
                GO TO 15
             ENDIF
         ENDIF
      ELSE
         WRITE(CHMAIL,4000) IEVENT
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GPHEAD: IDRUN=',I4,' IDEVT=',I7,
     &       ' Random Seeds = ',I10,2X,I10)
 1001 FORMAT(20X,10I10)
 2000 FORMAT(14X,'User words')
 3000 FORMAT( 1X,9(2X,E12.5))
 4000 FORMAT(' **** GPHEAD: No Header bank in memory for ievent:',I7)
*
      END

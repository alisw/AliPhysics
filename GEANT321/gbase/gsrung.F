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
      SUBROUTINE GSRUNG (NWBUF,UBUF,IADR)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Stores User's informations in RUNG user bank.            *
C     *       IADR is the stating location in the bank                 *
C     *                                                                *
C     *       on output IADR is the new starting location              *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author     M.Maire                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
*
      DIMENSION UBUF(1)
*
      IF (JRUNG.GT.0) THEN
*
         NDATA = IADR + NWBUF
         IF (LQ(JRUNG-1).EQ.0) THEN
            CALL MZBOOK (IXCONS,JRU,JRUNG,-1,'RUNU',0,0,NDATA,3,0)
         ELSE
            JRU = LQ(JRUNG-1)
            NFREE = IQ(JRU-1) - NDATA
            IF (NFREE.LT.0) CALL MZPUSH (IXCONS,JRU,0,-NFREE,'I')
         ENDIF
*
         DO 10 I=1,NWBUF
            Q(JRU+IADR+I) = UBUF(I)
   10    CONTINUE
*
         IADR = IADR + NWBUF
*
      ELSE
         IADR = 0
         WRITE(CHMAIL,1000)
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GSRUNG: No RUNG bank existing in memory ')
*
      END

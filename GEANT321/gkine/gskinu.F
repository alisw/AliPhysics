*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GSKINU (ITRA,NWBUF,UBUF,IADR)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Stores User's informations in JKINE  user bank,          *
C     *       for track ITRA. (ITRA must be an already existing track) *
C     *       IADR is the stating location within the bank.            *
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
#include "geant321/gcflag.inc"
*
      DIMENSION UBUF(1)
*
      IF ((JKINE.GT.0).AND.(LQ(JKINE-ITRA).GT.0)) THEN
*
         JK = LQ(JKINE-ITRA)
         NDATA = IADR + NWBUF
         IF (LQ(JK-1).EQ.0) THEN
            CALL MZBOOK (IXDIV,JKU,JK,-1,'KINU',0,0,NDATA,3,0)
            IQ(JKU-5) = ITRA
         ELSE
            JKU = LQ(JK-1)
            NFREE = IQ(JKU-1) - NDATA
            IF (NFREE.LT.0) CALL MZPUSH (IXDIV,JKU,0,-NFREE,'I')
         ENDIF
*
         DO 10 I=1,NWBUF
            Q(JKU+IADR+I) = UBUF(I)
   10    CONTINUE
*
         IADR = IADR + NWBUF
      ELSE
*
         IADR = 0
         WRITE(CHMAIL,1000) IEVENT,ITRA
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GSKINU: No JKINE bank in memory for ievent:',I7,
     &       '   or track:',I4,'does not exist')
*
      END

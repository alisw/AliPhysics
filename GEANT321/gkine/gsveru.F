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
      SUBROUTINE GSVERU (IVTX,NWBUF,UBUF,IADR)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Stores User's informations in VERTEX user bank,          *
C     *       for vertex IVTX. (IVTX must be an already existing vertex)*
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
      IF ((JVERTX.GT.0).AND.(LQ(JVERTX-IVTX).GT.0)) THEN
*
         JV = LQ(JVERTX-IVTX)
         NDATA = IADR + NWBUF
         IF (LQ(JV-1).EQ.0) THEN
            CALL MZBOOK (IXDIV,JVU,JV,-1,'VERU',0,0,NDATA,3,0)
            IQ(JVU-5) = IVTX
         ELSE
            JVU = LQ(JV-1)
            NFREE = IQ(JVU-1) - NDATA
            IF (NFREE.LT.0) CALL MZPUSH (IXDIV,JVU,0,-NFREE,'I')
         ENDIF
*
         DO 10 I=1,NWBUF
            Q(JVU+IADR+I) = UBUF(I)
   10    CONTINUE
*
         IADR = IADR + NWBUF
      ELSE
*
         IADR = 0
         WRITE(CHMAIL,1000) IEVENT,IVTX
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GSVERU: No Vertex bank in memory for ievent:',I7,
     &       '   or vertex:',I4,'does not exist')
*
      END

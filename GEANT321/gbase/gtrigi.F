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
      SUBROUTINE GTRIGI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Initialises event partition                              *
C.    *    ==>Called by : GRUN ,<GXINT> GINC4                          *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvdma.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctrak.inc"
C.
C.    ------------------------------------------------------------------
C.
      IEOTRI=0
      NTRACK=0
      NVERTX=0
      IDEBUG=0
      TOFG  =0.
C
C               Print event number and random number generator
C
      IF(RAYTRA.NE.1.)CALL GRNDMQ(NRNDM(1),NRNDM(2),0,'G')
C
C               Create event header bank
C
      CALL MZBOOK(IXDIV,JHEAD,JHEAD,1,'HEAD', 1, 1,NHEAD,2,0)
      IDEVT=IDEVT+1
      IQ(JHEAD+1)=IDRUN
      IQ(JHEAD+2)=IDEVT
      IQ(JHEAD+3)=NRNDM(1)
      IQ(JHEAD+4)=NRNDM(2)
C
      IF(ITEST.GT.0)THEN
         IF(MOD(IEVENT,ITEST).EQ.0)THEN
           IF(RAYTRA.NE.1.)THEN
            WRITE (CHMAIL,1000) IEVENT,IDEVT,(NRNDM(I),I = 1,2)
            CALL GMAIL(0,0)
           ENDIF
         ENDIF
      ENDIF
C
      IF(IEVENT.GE.IDEMIN.AND.IEVENT.LE.IDEMAX)IDEBUG=1
C
 1000 FORMAT(' **** GTRIGI: IEVENT=',I7,' IDEVT=',I7,
     +' Random Seeds = ',I10,2X,I10)
  99  RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GGET (LUN,KEYSU,NUKEYS,IDENT,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to read in data structures                       *
C.    *                                                                *
C.    *       LUN      Logical unit number                             *
C.    *       KEYSU    Keywords to select data structures              *
C.    *       NKEYS    Number of keywords                              *
C.    *       IER      Error flag                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT,GUKINE                        *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gccuts.inc"
      COMMON/QUEST/IQUEST(100)
C
      CHARACTER*4 KLEY(22)
      CHARACTER*4 KEYSU(1)
      DIMENSION KEYS(22),IUHEAD(2)
      DIMENSION KSEL(14),LKEY(22),LKNUM(22),LINK(14),JLINK(17)
      EQUIVALENCE (JLINK(1),JDIGI)
      SAVE IFIRST,LKEY
      DATA LINK/7,6,13,16,8,10,2,9,3,15,5,17,4,1/
      DATA KLEY/'PART','MATE','TMED','VOLU','ROTM','SETS','DRAW','RUNG'
     +         ,'INIT','INIT','INIT','INIT','INIT','INIT','INIT','INIT'
     +         ,'HEAD','KINE','KINE','JXYZ','HITS','DIGI'/
      DATA LKNUM/1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,13,14/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH(KLEY,LKEY,4,88)
      ENDIF
*
      WRITE(CHMAIL,10000)
      CALL GMAIL(0,0)
10000 FORMAT(' *** GGET *** Obsolete routine. Please use GFIN')
C
      IDENT=-1
      IER    = 0
      NKEYS=IABS(NUKEYS)
      IF (NKEYS.LE.0)                                 GO TO 99
      CALL UCTOH(KEYSU,KEYS,4,4*NKEYS)
C
      IF(NUKEYS.LT.0)THEN
        I1=1
        I2=15
        K1=1
        K2=7
      ELSE
        I1=18
        I2=22
        K1=10
        K2=14
      ENDIF
C
      NKT=0
      DO 10 K=K1,K2
  10  KSEL(K)=0
      DO 20 I=I1,I2
      N=LKNUM(I)
      DO 20 IK=1,NKEYS
      IF(KEYS(IK).EQ.LKEY(I))THEN
         KSEL(N)=1
         NKT=NKT+1
      ENDIF
  20  CONTINUE
      IF(NKT.EQ.0)GO TO 99
      NUH=2
C
C               Go for next start of event data structure
C
      IF(NUKEYS.LT.0)THEN
         IF(JRUNG.NE.0)CALL MZDROP(IXCONS,JRUNG,' ')
         CALL FZIN(LUN,IXCONS,JRUNG,1,'E',NUH,IUHEAD)
         IF(IQUEST(1).GT.2)GO TO 90
         IDIV=IXCONS
      ELSE
         IF(JHEAD.NE.0)CALL MZDROP(IXDIV,JHEAD,' ')
         CALL FZIN(LUN,IXDIV,JHEAD,1,'E',NUH,IUHEAD)
         IF(IQUEST(1).GT.2)GO TO 90
         IDIV=IXDIV
      ENDIF
C
      IDENT= IUHEAD(1)
      NK   = IUHEAD(2)
      IF(NK.LE.0)GO TO 99
      IF(NK.GT.10)GO TO 99
      DO 30 I=1,NK
C
C              Read next header
C
         NUH=2
         CALL FZIN(LUN,IDIV,0,0,'S',NUH,IUHEAD)
         IF(IQUEST(1).GT.2)GO TO 90
         KS=IUHEAD(1)
         IF(KS.LE.0)GO TO 30
         IF(KS.GT.14)GO TO 30
         IF(KSEL(KS).EQ.0)GO TO 30
         IL=LINK(KS)
         IF(JLINK(IL).NE.0)CALL MZDROP(IDIV,JLINK(IL),' ')
C
C              Read pending data structure
C
         CALL FZIN(LUN,IDIV,JLINK(IL),1,'A',NUH,IUHEAD)
         IF(IQUEST(1).GT.2)GO TO 90
  30  CONTINUE
C
C             Fill header bank
C             Reconstruct NKVIEW,NVOLUM,NVERTX,NTRACK
C             Reconstruct NMATE, NTMED, NPART
C
      IF(NUKEYS.LT.0)THEN
         IF(KSEL(1).NE.0.AND.JPART.GT.0) NPART=IQ(JPART-2)
         IF(KSEL(2).NE.0.AND.JMATE.NE.0) NMATE=IQ(JMATE-2)
         IF(KSEL(3).NE.0.AND.JTMED.NE.0) THEN
            CALL UCOPY(Q(JTMED+1),CUTGAM,10)
            NTMED=IQ(JTMED-2)
         ENDIF
         IF(KSEL(4).NE.0.AND.JVOLUM.GT.0) THEN
            NVOLUM=0
            DO 40 J=1, IQ(JVOLUM-2)
               IF(LQ(JVOLUM-J).EQ.0) GO TO 50
               NVOLUM=NVOLUM+1
  40        CONTINUE
  50        CONTINUE
         END IF
         IF(KSEL(7).NE.0.AND.JDRAW.GT.0) NKVIEW=IQ(JDRAW-2)
      ENDIF
C
      IF(JHEAD.GT.0)THEN
         IDRUN=IQ(JHEAD+1)
         IDEVT=IQ(JHEAD+2)
      ENDIF
C
      IF(KSEL(10).GT.0)THEN
         NVERTX=0
         NTRACK=0
         IF(JVERTX.GT.0)NVERTX=IQ(JVERTX+1)
         IF(JKINE .GT.0)NTRACK=IQ(JKINE +1)
      ENDIF
      GO TO 99
C
C             Error, EOF,etc
C
  90  IER=IQUEST(1)
C
  99  RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GSAVE(LUN,KEYSU,NUKEYS,IDENT,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to write out data structures                     *
C.    *                                                                *
C.    *       LUN      Logical unit number                             *
C.    *       KEYSU    Keywords to select data structures              *
C.    *       NKEYS    Number of keywords                              *
C.    *       IER      Error flag                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT,GUOUT                         *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 KLEY(22)
      CHARACTER*4 KEYSU(1)
      DIMENSION KEYS(22),IUHEAD(2)
      DIMENSION KSEL(14),LKEY(22),LKNUM(22),LINK(14),JLINK(17)
      EQUIVALENCE (JLINK(1),JDIGI)
      COMMON/QUEST/IQUEST(100)
      SAVE IFIRST,LKEY
C
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
10000 FORMAT(' *** GSAVE *** Obsolete routine. Please use GFOUT')
C
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
      DO 10 K=K1,K2
  10  KSEL(K)=0
      NK=0
      DO 25 I=I1,I2
         N=LKNUM(I)
         DO 20 IK=1,NKEYS
            IF(KEYS(IK).EQ.LKEY(I))THEN
               IL=LINK(N)
               IF(JLINK(IL).NE.0)THEN
                  KSEL(N)=1
                  NK=NK+1
               ENDIF
            ENDIF
  20     CONTINUE
  25  CONTINUE
C
      IUHEAD(1)=IDENT
      IUHEAD(2)=NK
      IF(NUKEYS.LT.0)THEN
C
C======>      Write RUN header and constants
C
         CALL FZOUT(LUN,IXCONS,JRUNG,1,'L',2,2,IUHEAD)
         IF(IQUEST(1).NE.0)GO TO 90
         DO 30 I=1,7
            IF(KSEL(I).NE.0)THEN
               IL=LINK(I)
               CALL FZOUT(LUN,IXCONS,JLINK(IL),0,'L',2,1,I)
               IF(IQUEST(1).NE.0)GO TO 90
            ENDIF
  30     CONTINUE
C
      ELSE
C
C======>      Write event header and data structures
C             Released unused space in JHITS and JDIGI
C
         IF(KSEL(13).NE.0)CALL GRLEAS(JHITS)
         IF(KSEL(14).NE.0)CALL GRLEAS(JDIGI)
C
         CALL FZOUT(LUN,IXDIV,JHEAD,1,' ',2,2,IUHEAD)
         IF(IQUEST(1).NE.0)GO TO 90
         DO 40 I=10,14
            IF(KSEL(I).NE.0)THEN
               IL=LINK(I)
               CALL FZOUT(LUN,IXDIV ,JLINK(IL),0,'L',2,1,I)
               IF(IQUEST(1).NE.0)GO TO 90
            ENDIF
  40     CONTINUE
C
      ENDIF
      GO TO 99
C
C             Error
C
  90  IER=IQUEST(1)
C
  99  RETURN
      END

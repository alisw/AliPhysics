*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GRGET (KEYSU,ID1,ID2,ID3,ICYCLE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to read  data structures from a RZ file          *
C.    *                                                                *
C.    *       KEYSU    Keyword to select data structure                *
C.    *       ID1      First  RZ KEY identifier (ex IDRUN)             *
C.    *       ID2      Second RZ KEY identifier (ex IDEVT)             *
C.    *       ID3      Third  RZ KEY identifier (user free)            *
C.    *       ICYCLE   Cycle number                                    *
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
C
      DIMENSION CUTVEC(10)
      EQUIVALENCE (CUTVEC,CUTGAM)
C
      CHARACTER*4 KLEY(19)
      CHARACTER*4 KEYSU(1)
      DIMENSION KEYRZ(4)
      DIMENSION LINIT(8),LKINE(2),LTRIG(6)
      DIMENSION LKEY(19),LINK(10),JLINK(17),KSEL(17),IXD(17)
      EQUIVALENCE (JLINK(1),JDIGI)
      EQUIVALENCE (LKEY(18),KINIT),(LKEY(19),KTRIG),(LKEY(5),KKINE)
C
      SAVE IFIRST,LKEY
      DATA KLEY/'DIGI','DRAW','HEAD','HITS','KINE','MATE','PART','ROTM'
     +         ,'RUNG','SETS','STAK','STAT','TMED','TRAC','VERT','VOLU'
     +         ,'JXYZ','INIT','TRIG'/
      DATA IXD/2,1,2,2,2,8*1,2,2,1,2/
      DATA LINIT/2,6,7,8,9,10,13,16/
      DATA LKINE/5,15/
      DATA LTRIG/1,3,4,5,15,17/
      DATA NLINIT/8/
      DATA NLKINE/2/
      DATA NLTRIG/6/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH(KLEY,LKEY,4,76)
      ENDIF
*
      WRITE(CHMAIL,10000)
      CALL GMAIL(0,0)
10000 FORMAT(' *** GRGET *** Obsolete routine. Please use GRIN')
C
      CALL UCTOH(KEYSU,KEY,4,4)
      CALL VZERO(KSEL,17)
C
      IF(KEY.EQ.KINIT)THEN
         CALL UCOPY(LINIT,LINK,NLINIT)
         NLINK=NLINIT
      ELSEIF(KEY.EQ.KKINE)THEN
         CALL UCOPY(LKINE,LINK,NLKINE)
         NLINK=NLKINE
      ELSEIF(KEY.EQ.KTRIG)THEN
         CALL UCOPY(LTRIG,LINK,NLTRIG)
         NLINK=NLTRIG
      ELSE
         IL=IUCOMP(KEY,LKEY,17)
         IF(IL.EQ.0)GO TO 50
         LINK(1)=IL
         NLINK=1
      ENDIF
C
C               Read data structure(s)
C
      KEYRZ(2)=ID1
      KEYRZ(3)=ID2
      KEYRZ(4)=ID3
C
      DO 10 I=1,NLINK
         IL=LINK(I)
         IF(JLINK(IL).NE.0)THEN
            CALL MZDROP(IXCONS,JLINK(IL),' ')
            JLINK(IL)=0
         ENDIF
         KEYRZ(1)=LKEY(IL)
         IF(IXD(IL).NE.2)THEN
            IDIV=IXCONS
         ELSE
            IDIV=IXDIV
         ENDIF
         CALL RZIN(IDIV,JLINK(IL),1,KEYRZ,ICYCLE,' ')
         KSEL(IL)=1
   10 CONTINUE
C
C             Fill header bank
C             Reconstruct NKVIEW,NVOLUM,NVERTX,NTRACK
C             Reconstruct NMATE, NTMED, NPART
C
      IF(KSEL(6).NE.0.AND.JMATE.GT.0) NMATE=IQ(JMATE-2)
      IF(KSEL(7).NE.0.AND.JPART.GT.0) NPART=IQ(JPART-2)
      IF(KSEL(13).NE.0.AND.JTMED.NE.0 ) THEN
         DO 20 J=1,10
            CUTVEC(J) = Q(JTMED+J)
   20    CONTINUE
         NTMED=IQ(JTMED-2)
      ENDIF
      IF(KSEL(16).NE.0.AND.JVOLUM.GT.0) THEN
         NVOLUM=0
         DO 30 J=1, IQ(JVOLUM-2)
            IF(LQ(JVOLUM-J).EQ.0) GO TO 40
            NVOLUM=NVOLUM+1
   30    CONTINUE
   40    CONTINUE
      END IF
      IF(KSEL( 2).NE.0.AND.JDRAW.GT.0 ) NKVIEW=IQ(JDRAW-2)
C
      IF(JHEAD.GT.0)THEN
         IDRUN=IQ(JHEAD+1)
         IDEVT=IQ(JHEAD+2)
      ENDIF
C
      IF(KSEL(5).GT.0)THEN
         NVERTX=0
         NTRACK=0
         IF(JVERTX.GT.0)NVERTX=IQ(JVERTX+1)
         IF(JKINE .GT.0)NTRACK=IQ(JKINE +1)
      ENDIF
C
   50 RETURN
      END

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
      SUBROUTINE GRSAVE(KEYSU,ID1,ID2,ID3,ICYCLE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to write out data structures on a RZ file        *
C.    *                                                                *
C.    *       KEYSU    Keyword to select data structure(s)             *
C.    *       ID1      First  RZ KEY identifier (ex IDRUN)             *
C.    *       ID2      Second RZ KEY identifier (ex IDEVT)             *
C.    *       ID3      Third  RZ KEY identifier (user free)            *
C.    *       ICYCLE   Cycle number (output)                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT,GUOUT                         *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 KLEY(19)
      CHARACTER*4 KEYSU(1)
      DIMENSION KEYRZ(4)
      DIMENSION LINIT(8),LKINE(2),LTRIG(6)
      DIMENSION LKEY(19),LINK(10),JLINK(17)
      EQUIVALENCE (JLINK(1),JDIGI)
      EQUIVALENCE (LKEY(18),KINIT),(LKEY(19),KTRIG),(LKEY(5),KKINE)
      SAVE IFIRST,LKEY
C
      DATA KLEY/'DIGI','DRAW','HEAD','HITS','KINE','MATE','PART','ROTM'
     +         ,'RUNG','SETS','STAK','STAT','TMED','TRAC','VERT','VOLU'
     +         ,'JXYZ','INIT','TRIG'/
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
10000 FORMAT(' *** GRSAVE *** Obsolete routine. Please use GROUT')
C
      CALL UCTOH(KEYSU,KEY,4,4)
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
         IF(IL.EQ.0)GO TO 99
         LINK(1)=IL
         NLINK=1
      ENDIF
C
C               Write data structure(s)
C
      KEYRZ(2)=ID1
      KEYRZ(3)=ID2
      KEYRZ(4)=ID3
C
      DO 10 I=1,NLINK
         IL=LINK(I)
         IF(JLINK(IL).EQ.0)GO TO 10
         KEYRZ(1)=LKEY(IL)
         IF(IL.EQ.4)CALL GRLEAS(JHITS)
         IF(IL.EQ.1)CALL GRLEAS(JDIGI)
         CALL RZOUT(IXSTOR,JLINK(IL),KEYRZ,ICYCLE,'L')
  10  CONTINUE
C
  99  RETURN
      END

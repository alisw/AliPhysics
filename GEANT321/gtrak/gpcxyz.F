*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 06/12/94  16.50.58  by  S.Ravndal
*-- Author :
      SUBROUTINE GPCXYZ
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        Print track and volume parameters at current point      *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      DIMENSION MECNAM(20)
      CHARACTER*4 KUN1,KUN2
      SAVE IEVOLD,NTMOLD
      DATA IEVOLD,NTMOLD/-1,-1/
C.
C.    ------------------------------------------------------------------
C.
      IF(IFINIT(9).EQ.0)THEN
         IFINIT(9)=1
         IEVOLD=-1
         NTMOLD=-1
      ENDIF
C
      NM=NMEC
      IF(NM.EQ.0)THEN
         MECNAM(1)=NAMEC(29)
         NM=1
      ELSE
         DO 10 I=1,NMEC
            MEC=LMEC(I)
            IF(MEC.LE.MAXMEC) THEN
               MECNAM(I)=NAMEC(MEC)
            ELSEIF(MEC-100.LE.MAXME1.AND.MEC-100.GT.0) THEN
               MECNAM(I)=NAMEC1(MEC-100)
            ENDIF
  10     CONTINUE
      ENDIF
C
      IF(IEVENT.EQ.IEVOLD.AND.NTMULT.EQ.NTMOLD)GO TO 20
C
C
      TOFGN=TOFG*1.E+9
      WRITE(CHMAIL,1000)ITRA,ISTAK,NTMULT,(NAPART(I),I=1,5),TOFGN
      CALL GMAIL(0,0)
      WRITE(CHMAIL,1100)
      CALL GMAIL(0,0)
      IEVOLD=IEVENT
      NTMOLD=NTMULT
C
  20  R=SQRT(VECT(1)**2+VECT(2)**2)
      CALL GEVKEV(DESTEP,DESU,KUN1)
      CALL GEVKEV(GEKIN ,GEKU,KUN2)
      NS = 1
      NW = MIN(NS+5,NM)
      WRITE(CHMAIL,2000)(VECT(I),I=1,3),R,NAMES(NLEVEL),NUMBER(NLEVEL)
     +      ,SLENG,STEP,DESU,KUN1,GEKU,KUN2,(MECNAM(I),I=NS,NW)
  30  CALL GMAIL(0,0)
      IF(NM.GT.NW) THEN
         NS = NW + 1
         NW = MIN(NS+5,NM)
         WRITE(CHMAIL,3000) (MECNAM(I),I=NS,NW)
         GOTO 30
      ENDIF
C
 1000 FORMAT(' =====> TRACK ',I4,' STACK NR',I4,' NTMULT=',I5,5X,
     +5A4,5X,'TOFG =',F10.3,' NS')
 1100 FORMAT('       X          Y          Z          R   ',
     +'  NAME  NUMBER',
     +'   SLENG      STEP      DESTEP     GEKIN   MECHANISMS')
 2000 FORMAT(1X,4F11.4,1X,A4,1X,I4,1X,2F10.4,F7.1,A4,F9.3,A4,1X,
     +          6(A4,1X))
 3000 FORMAT(101X,6(A4,1X))
C
      END

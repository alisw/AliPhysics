#include "isajet/pilot.h"
      SUBROUTINE ORDER(ID,MODEIN,MODOUT,MEOUT)
C
C          Search for mode MODEIN of particle ID in standard /DKYTAB/.
C          If found, return MODOUT = standard order and MEOUT=MELEM.
C          Otherwise return MODOUT = MODEIN and MEOUT=0.
C          If ID<0, use antiparticles instead.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/itapes.inc"
#include "isajet/dkytab.inc"
#include "isajet/force.inc"
C
      INTEGER ID,MODEIN(5),MODOUT(5),MODTST(5)
      INTEGER IFL1,IFL2,IFL3,JSPIN,INDEX,LOOK0,IUSE(5),ISAME,I,J,
     $NADD,NADDI,K,K1,K2,IDANTI,MEOUT
C
C          Find standard starting point
C
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IF(LOOK(INDEX).GT.0) THEN
        LOOK0=LOOK(INDEX)
      ELSEIF(LOOK(INDEX).LT.0) THEN
        LOOK0=LOOKST(-LOOK(INDEX))
      ELSE
        ISAME=0
        GO TO 300
      ENDIF
C
C          Find NADD
C
      DO 100 I=1,5
100   IF(MODEIN(I).NE.0) NADD=I
C
C          If ID<0, compare antiparticles
C
      IF(ID.GE.0) THEN
        DO 110 K=1,NADD
110     MODTST(K)=MODEIN(K)
      ELSE
        DO 120 K=1,NADD
120     MODTST(K)=IDANTI(MODEIN(K))
      ENDIF
C
C          Scan all modes starting at LOOK0. Check for correct NADD.
C          Then check that particles match in arbitrary order.
C
      IF(LOOK0.LE.0) GO TO 300
      DO 200 I=LOOK0,MXDKY
        DO 210 K=1,5
210     IF(MODE(K,I).NE.0) NADDI=K
        IF(NADDI.EQ.NADD) THEN
          DO 220 K=1,5
220       IUSE(K)=0
C
          DO 230 K1=1,NADD
            DO 240 K2=1,NADD
              IF(MODTST(K1).EQ.MODE(K2,I).AND.IUSE(K2).EQ.0) THEN
                IUSE(K2)=K1
                GO TO 230
              ENDIF
240         CONTINUE
            GO TO 201
230       CONTINUE
C
          ISAME=I
          GO TO 300
        ENDIF
201     IF(CBR(I).GE.1.) THEN
          ISAME=0
          GO TO 300
        ENDIF
200   CONTINUE
      STOP 99
C
C          Return matching mode or original mode.
C
300   IF(ISAME.EQ.0) THEN
        WRITE(ITLIS,3001)
3001    FORMAT(' ***** WARNING: NONSTANDARD MODE')
        DO 310 K=1,5
310     MODOUT(K)=MODEIN(K)
        MEOUT=0
      ELSEIF(ID.GT.0) THEN
        DO 320 K=1,5
320     MODOUT(K)=MODE(K,ISAME)
        MEOUT=MELEM(ISAME)
      ELSE
        DO 330 K=1,5
330     MODOUT(K)=IDANTI(MODE(K,ISAME))
        MEOUT=MELEM(ISAME)
      ENDIF
C
      RETURN
      END

#include "isajet/pilot.h"
      SUBROUTINE PRTLST(JTLIS,AMY,AMX)
C
C          List defined particles. AMY, AMX are the masses of the 
C          fourth generation quarks. If a negative mass is given,
C          then these are not listed.
C          This must be linked with ISAJET, including ALDATA.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/qlmass.inc"
#include "isajet/wcon.inc"
      INTEGER JTLIS
      REAL AMY,AMX
      INTEGER IFL1,IFL2,IFL3,JSPIN,INDEX,I,ID
      REAL AM,CG,AMASS,CHARGE
      CHARACTER*8 LB,LABEL
C
C          Initialize SUSY masses to 0. Remember offset of 1 from KL.
C
      DO 100 I=22,NQLEP
        AMLEP(I)=0.
100   CONTINUE
      AMLEP(7)=AMY
      AMLEP(8)=AMX
      CALL FLAVOR(80,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=WMASS(2)
      CALL FLAVOR(90,IFL1,IFL2,IFL3,JSPIN,INDEX)
      AMLEP(INDEX)=WMASS(4)
      WRITE(JTLIS,101) AMY,AMX
C
C          Loop over IDENT's
C
      DO 200 I=1,40000
        ID = I
        CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
        IF(AMX.LT.0..OR.AMY.LT.0.) THEN
          IF(IABS(IFL1).GT.6.OR.IABS(IFL2).GT.6.OR.IABS(IFL3).GT.6)
     $    GO TO 200
        ENDIF
        IF(INDEX.GT.0) THEN
          LB = LABEL(ID)
          IF(LB.NE.'ERR') THEN
            AM = AMASS(ID)
            CG = CHARGE(ID)
            WRITE(JTLIS,102) ID,LB,AM,CG,IFL1,IFL2,IFL3,JSPIN,INDEX
          ENDIF
        ENDIF
        ID = -I
        CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
C          Eliminate bad ID's:
        IF(INDEX.GT.0) THEN
          LB = LABEL(ID)
          IF(LB.NE.'ERR') THEN
            AM = AMASS(ID)
            CG = CHARGE(ID)
            WRITE(JTLIS,102) ID,LB,AM,CG,IFL1,IFL2,IFL3,JSPIN,INDEX
          ENDIF
        ENDIF
 200  CONTINUE
      RETURN
 101  FORMAT(10X,'ISAJET PARTICLES, M(Y) =',F10.3,'  M(X) =',F10.3//
     1 5X,'ID',4X,'PARTICLE',8X,'MASS',4X,'CHARGE',
     2 4X,'---FLAVOR---',4X,'SPIN',4X,'INDEX')
 102  FORMAT(1X,I6,4X,A8,F12.6,F10.2,4X,3I4,I8,I9)
      END

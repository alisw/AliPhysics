#include "isajet/pilot.h"
      LOGICAL FUNCTION SETTYP(LPRT)
C
C          Set JETTYPE flags and WMODE flags for WPAIR.
C          Set WMODES and ZMODES flags for secondary W+- and Z0.
C          Return .FALSE. if no error, .TRUE. otherwise.
C
C          Ver 7.18: Initialize all GOQ to false (limit = MXGOQ)
C                    Use LISTSS for Higgs if GOMSSM
C          Ver 7.29: SUSY Higgs decays are done in SETHSS and SIGHSS
C                    using LISTSS order, so SUSY list should be used.
C                    I.e., 7.18 fix was wrong.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/keys.inc"
#include "isajet/types.inc"
#include "isajet/q1q2.inc"
#include "isajet/xmssm.inc"
C
      INTEGER JET,K,I,IW,LPRT
      INTEGER NLIST
      CHARACTER*8 WORD,BLANK,LIST(30),LISTW(4),LISTXY(4),LISTSS(85)
      DATA BLANK/'        '/
      DATA LIST/'GL','UP','UB','DN','DB','ST','SB','CH','CB','BT','BB',
     $'TP','TB','NUE','ANUE','E-','E+','NUM','ANUM','MU-','MU+',
     $'NUT','ANUT','TAU-','TAU+','GM','W+','W-','Z0','HIGGS'/
      DATA LISTW/'GM','W+','W-','Z0'/
      DATA LISTXY/'Y','YB','X','XB'/
      DATA LISTSS/'GLSS',
     $'UPSSL','UBSSL','DNSSL','DBSSL','STSSL','SBSSL','CHSSL','CBSSL',
     $'BTSS1','BBSS1','TPSS1','TBSS1',
     $'UPSSR','UBSSR','DNSSR','DBSSR','STSSR','SBSSR','CHSSR','CBSSR',
     $'BTSS2','BBSS2','TPSS2','TBSS2',
     $'W1SS+','W1SS-','W2SS+','W2SS-','Z1SS','Z2SS','Z3SS','Z4SS',
     $'NUEL','ANUEL','EL-','EL+','NUML','ANUML','MUL-','MUL+',
     $'NUTL','ANUTL','TAU1-','TAU1+','ER-','ER+','MUR-','MUR+',
     $'TAU2-','TAU2+',
     $'GL','UP','UB','DN','DB','ST','SB','CH','CB','BT','BB',
     $'TP','TB','NUE','ANUE','E-','E+','NUM','ANUM','MU-','MU+',
     $'NUT','ANUT','TAU-','TAU+','GM','W+','W-','Z0',
     $'HL0','HH0','HA0','H+','H-'/
C
      SETTYP=.FALSE.
C       
      IF(KEYS(5)) GO TO 5
      IF(KEYS(2).AND.GOMSSM) GO TO 5
      IF(KEYS(6).OR.KEYS(9)) GO TO 6
      IF(KEYS(7).AND..NOT.GOMSSM) GO TO 7
      IF(KEYS(7).AND.GOMSSM) GO TO 5
      IF(KEYS(10).AND.GOMSSM) GO TO 5
C
C          JETTYPE flags all processes except WPAIR and HIGGS.
C          NJTTYP is set in READIN to number of non-blank values read.
C          Check for legal jet type names and set appropriate flags.
C
      DO 1000 JET=1,MXGOJ
        IF(NJTTYP(JET).EQ.0) GO TO 1000
C          Initialize everything to .FALSE.
        GOALL(JET)=.FALSE.
        DO 1100 K=1,MXGOQ
          GOQ(K,JET)=.FALSE.
1100    CONTINUE
C          Loop over non-blank JETTYPE entries
        DO 1200 I=1,NJTTYP(JET)
          WORD=JETYP(I,JET)
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 1200
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            GOALL(JET)=.TRUE.
            DO 1210 K=1,MXGOQ
1210        GOQ(K,JET)=.TRUE.
            GO TO 1000
          ENDIF
C          Quarks
          IF(WORD.EQ.'QUARKS  ') THEN
            DO 1220 K=2,13
1220        GOQ(K,JET)=.TRUE.
            GO TO 1200
          ENDIF
C          Charged leptons
          IF(WORD.EQ.'LEPTONS ') THEN
            DO 1230 K=16,24,4
            GOQ(K,JET)=.TRUE.
1230        GOQ(K+1,JET)=.TRUE.
            GO TO 1200
          ENDIF
C          Neutrinos
          IF(WORD.EQ.'NUS     ') THEN
            DO 1240 K=14,22,4
            GOQ(K,JET)=.TRUE.
1240        GOQ(K+1,JET)=.TRUE.
            GO TO 1200
          ENDIF
C          Explicit types
C          E+E- now also contains W+, W-, Z0
          IF(KEYS(10).OR.KEYS(11).OR.KEYS(12)) THEN
            NLIST=30
          ELSE IF(KEYS(2)) THEN
            NLIST=29
          ELSE
            NLIST=25
          ENDIF
          DO 1250 K=1,NLIST
            IF(WORD.EQ.LIST(K)) THEN
              GOQ(K,JET)=.TRUE.
              GO TO 1200
            ENDIF
1250      CONTINUE
C          Special types for TWOJET
          DO 1270 K=1,4
            IF(KEYS(1).AND.WORD.EQ.LISTXY(K)) THEN
              GOQ(13+K,JET)=.TRUE.
              GO TO 1200
            ENDIF
1270      CONTINUE
C          Special type for PHOTON
          IF(KEYS(8).AND.WORD.EQ.LISTW(1)) THEN
            GOQ(26,JET)=.TRUE.
            GO TO 1200
          ENDIF
C          Error
          WRITE(ITLIS,1300) WORD,JET
1300      FORMAT(1X,A8,' IS NOT RECOGNIZABLE FOR JETTYPE',I1)
          SETTYP=.TRUE.
1200    CONTINUE
1000  CONTINUE
      GO TO 4000
C
C          JETTYPE flags for SUSY
C
5     DO 5000 JET=1,2
        IF(NJTTYP(JET).EQ.0) GO TO 5000
        GOALL(JET)=.FALSE.
        DO 5100 K=1,MXGOQ
5100    GOQ(K,JET)=.FALSE.
        DO 5200 I=1,NJTTYP(JET)
          WORD=JETYP(I,JET)
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 5200
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            GOALL(JET)=.TRUE.
            DO 5210 K=1,85
5210        GOQ(K,JET)=.TRUE.
            GO TO 5000
          ENDIF
C          Squarks
          IF(WORD.EQ.'SQUARKS ') THEN
            DO 5220 K=2,25
5220        GOQ(K,JET)=.TRUE.
            GO TO 5200
          ENDIF
C           Gauginos
          IF(WORD.EQ.'GAUGINOS') THEN
            DO 5230 K=26,33
5230        GOQ(K,JET)=.TRUE.
            GO TO 5200
          ENDIF
C           Sleptons
          IF(WORD.EQ.'SLEPTONS') THEN
            DO 5240 K=34,51
5240        GOQ(K,JET)=.TRUE.
            GO TO 5200
          ENDIF
C          Explicit susy types
          DO 5300 K=1,85
            IF(WORD.EQ.LISTSS(K)) THEN
              GOQ(K,JET)=.TRUE.
              GO TO 5200
            ENDIF
5300      CONTINUE
5200    CONTINUE
5000  CONTINUE
      GO TO 4000
C
C          JETTYPE and WMODE flags for WPAIR
C          NJTTYP and NWWTYP are the number of non-blank values.
C
6     DO 2000 JET=1,2
        IF(NJTTYP(JET).EQ.0) GO TO 2300
C          Initialize to FALSE
        GOALL(JET)=.FALSE.
        DO 2100 K=1,4
2100    GOQ(K,JET)=.FALSE.
C
C          Loop over non-blank JETTYPE flags
C
        DO 2200 I=1,NJTTYP(JET)
          WORD=JETYP(I,JET)
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 2200
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            GOALL(JET)=.TRUE.
            DO 2210 K=1,4
2210        GOQ(K,JET)=.TRUE.
            GO TO 2300
          ENDIF
C          Explicit types
          DO 2220 K=1,4
            IF(WORD.EQ.LISTW(K)) THEN
              GOQ(K,JET)=.TRUE.
              GO TO 2200
            ENDIF
2220      CONTINUE
C          Error
          WRITE(ITLIS,1300) WORD,JET
          SETTYP=.TRUE.
2200    CONTINUE
C
C          Loop over nonblank WMODE flags
C
2300    IF(NWWTYP(JET).EQ.0) GO TO 2000
        ALLWW(JET)=.FALSE.
C         Initialize everything to FALSE
        DO 2350 K=1,25
2350    GOWW(K,JET)=.FALSE.
C
        DO 2400 I=1,NWWTYP(JET)
          WORD=WWTYP(I,JET)
          IF(WORD.NE.BLANK) NWWTYP(JET)=I
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 2400
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            ALLWW(JET)=.TRUE.
            DO 2410 K=1,25
2410        GOWW(K,JET)=.TRUE.
            GO TO 2000
          ENDIF
C          Quarks
          IF(WORD.EQ.'QUARKS  ') THEN
            DO 2420 K=2,13
2420        GOWW(K,JET)=.TRUE.
            GO TO 2400
          ENDIF
C          Charged leptons
          IF(WORD.EQ.'LEPTONS ') THEN
            DO 2430 K=16,24,4
            GOWW(K,JET)=.TRUE.
2430        GOWW(K+1,JET)=.TRUE.
            GO TO 2400
          ENDIF
C          Neutrinos
          IF(WORD.EQ.'NUS     ') THEN
            DO 2440 K=14,22,4
            GOWW(K,JET)=.TRUE.
2440        GOWW(K+1,JET)=.TRUE.
            GO TO 2400
          ENDIF
C          Explicit types
          DO 2450 K=1,25
            IF(WORD.EQ.LIST(K)) THEN
              GOWW(K,JET)=.TRUE.
              GO TO 2400
            ENDIF
2450      CONTINUE
C          Error
          WRITE(ITLIS,2500) WORD,JET
2500      FORMAT(1X,A8,' IS NOT A VALID CODE FOR WMODE',I1)
          SETTYP=.TRUE.
2400    CONTINUE
2000  CONTINUE
      GO TO 4000
C
C          JETTYPE and WMODE flags for HIGGS
C          SUSY HIGGS uses LISTSS order and hence SUSY part
C
7     DO 3000 JET=1,2
        IF(NJTTYP(JET).EQ.0) GO TO 3300
C          Initialize to FALSE
        GOALL(JET)=.FALSE.
        DO 3100 K=1,MXGOQ
3100    GOQ(K,JET)=.FALSE.
C
C          Loop over non-blank JETTYPE flags
C
        DO 3200 I=1,NJTTYP(JET)
          WORD=JETYP(I,JET)
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 3200
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            GOALL(JET)=.TRUE.
            DO 3210 K=1,MXGOQ
3210        GOQ(K,JET)=.TRUE.
            GO TO 3300
          ENDIF
C          Quarks
          IF(WORD.EQ.'QUARKS  ') THEN
            DO 3220 K=2,13
3220        GOQ(K,JET)=.TRUE.
            GO TO 3200
          ENDIF
C          Charged leptons
          IF(WORD.EQ.'LEPTONS ') THEN
            DO 3240 K=16,24,4
            GOQ(K,JET)=.TRUE.
3240        GOQ(K+1,JET)=.TRUE.
            GO TO 3200
          ENDIF
          DO 3250 K=1,85
            IF(WORD.EQ.LIST(K)) THEN
              GOQ(K,JET)=.TRUE.
              GO TO 3200
            ENDIF
3250      CONTINUE
C          Error
          WRITE(ITLIS,1300) WORD,JET
          SETTYP=.TRUE.
3200    CONTINUE
C
C          Loop over nonblank WMODE flags
C
3300    IF(NWWTYP(JET).EQ.0) GO TO 3000
        ALLWW(JET)=.FALSE.
C         Initialize everything to FALSE
        DO 3350 K=1,25
3350    GOWW(K,JET)=.FALSE.
C
        DO 3400 I=1,NWWTYP(JET)
          WORD=WWTYP(I,JET)
          IF(WORD.NE.BLANK) NWWTYP(JET)=I
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 3400
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            ALLWW(JET)=.TRUE.
            DO 3410 K=1,25
3410        GOWW(K,JET)=.TRUE.
            GO TO 3000
          ENDIF
C          Quarks
          IF(WORD.EQ.'QUARKS  ') THEN
            DO 3420 K=2,13
3420        GOWW(K,JET)=.TRUE.
            GO TO 3400
          ENDIF
C          Charged leptons
          IF(WORD.EQ.'LEPTONS ') THEN
            DO 3430 K=16,24,4
            GOWW(K,JET)=.TRUE.
3430        GOWW(K+1,JET)=.TRUE.
            GO TO 3400
          ENDIF
C          Neutrinos
          IF(WORD.EQ.'NUS     ') THEN
            DO 3440 K=14,22,4
            GOWW(K,JET)=.TRUE.
3440        GOWW(K+1,JET)=.TRUE.
            GO TO 3400
          ENDIF
C          Explicit types
          DO 3450 K=1,25
            IF(WORD.EQ.LIST(K)) THEN
              GOWW(K,JET)=.TRUE.
              GO TO 3400
            ENDIF
3450      CONTINUE
C          Error
          WRITE(ITLIS,2500) WORD,JET
3500      FORMAT(1X,A8,' IS NOT A VALID CODE FOR WMODE',I1)
          SETTYP=.TRUE.
3400    CONTINUE
3000  CONTINUE
C
C          Set WMODES and ZMODES flags for secondary W+- and Z0
C
4000  DO 4100 IW=1,3
        IF(NWMODE(IW).EQ.0) GO TO 4100
C           Initialize everything to .FALSE.
        DO 4200 K=1,25
4200    GOWMOD(K,IW)=.FALSE.
C          Loop over non-blank WMODE entries
        DO 4300 I=1,NWMODE(IW)
          WORD=WMODES(I,IW)
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 4300
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            DO 4310 K=1,25
4310        GOWMOD(K,IW)=.TRUE.
            GO TO 4100
          ENDIF
C          Quarks
          IF(WORD.EQ.'QUARKS  ') THEN
            DO 4320 K=2,13
4320        GOWMOD(K,IW)=.TRUE.
            GO TO 4300
          ENDIF
C          Charged leptons
          IF(WORD.EQ.'LEPTONS ') THEN
            DO 4330 K=16,24,4
            GOWMOD(K,IW)=.TRUE.
4330        GOWMOD(K+1,IW)=.TRUE.
            GO TO 4300
          ENDIF
C          Neutrinos
          IF(WORD.EQ.'NUS     ') THEN
            DO 4340 K=14,22,4
            GOWMOD(K,IW)=.TRUE.
4340        GOWMOD(K+1,IW)=.TRUE.
            GO TO 4300
          ENDIF
C          Explicit types
          DO 4350 K=1,25
            IF(WORD.EQ.LIST(K)) THEN
              GOWMOD(K,IW)=.TRUE.
              GO TO 4300
            ENDIF
4350      CONTINUE
C          Error
          WRITE(ITLIS,4380) WORD
4380      FORMAT(1X,A8,' IS NOT RECOGNIZABLE FOR SECONDARY WS')
          SETTYP=.TRUE.
4300    CONTINUE
4100  CONTINUE
C
C          Loop over nonblank WMODE flags FOR WHIGGS
C
      IF (KEYS(10)) THEN
      DO 6000 JET=1,2
6300    IF(NWWTYP(JET).EQ.0) GO TO 6000
        ALLWW(JET)=.FALSE.
C         Initialize everything to FALSE
        DO 6350 K=1,25
6350    GOWW(K,JET)=.FALSE.
C
        DO 6400 I=1,NWWTYP(JET)
          WORD=WWTYP(I,JET)
          IF(WORD.NE.BLANK) NWWTYP(JET)=I
C          Blank
          IF(WORD.EQ.BLANK) THEN
            GO TO 6400
          ENDIF
C          All
          IF(WORD.EQ.'ALL     ') THEN
            ALLWW(JET)=.TRUE.
            DO 6410 K=1,25
6410        GOWW(K,JET)=.TRUE.
            GO TO 6000
          ENDIF
C          Quarks
          IF(WORD.EQ.'QUARKS  ') THEN
            DO 6420 K=2,13
6420        GOWW(K,JET)=.TRUE.
            GO TO 6400
          ENDIF
C          Charged leptons
          IF(WORD.EQ.'LEPTONS ') THEN
            DO 6430 K=16,24,4
            GOWW(K,JET)=.TRUE.
6430        GOWW(K+1,JET)=.TRUE.
            GO TO 6400
          ENDIF
C          Neutrinos
          IF(WORD.EQ.'NUS     ') THEN
            DO 6440 K=14,22,4
            GOWW(K,JET)=.TRUE.
6440        GOWW(K+1,JET)=.TRUE.
            GO TO 6400
          ENDIF
C          Explicit types
          DO 6450 K=1,25
            IF(WORD.EQ.LIST(K)) THEN
              GOWW(K,JET)=.TRUE.
              GO TO 6400
            ENDIF
6450      CONTINUE
C          Error
          WRITE(ITLIS,6500) WORD,JET
6500      FORMAT(1X,A8,' IS NOT A VALID CODE FOR WMODE',I1)
          SETTYP=.TRUE.
6400    CONTINUE
6000  CONTINUE
      END IF
      RETURN
      END

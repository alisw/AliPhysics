#include "isajet/pilot.h"
      SUBROUTINE SETHSS
C
C          Set the MSSM Higgs parameters in /HCON/.
C          HMASS  = Higgs mass for HTYPE
C          HGAM   = Higgs width
C          HGAMSS = Higgs partial widths. Note HGAMSS is not
C                   necessarily diagonal for SUSY decays.
C          ZSTARS = minimum allowed mass for Z*
C
C          Note LISTSS(78) => W+, LISTSS(79) => W-, LISTSS(80) => Z0
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/hcon.inc"
#include "isajet/listss.inc"
#include "isajet/q1q2.inc"
#include "isajet/ssmode.inc"
#include "isajet/sstype.inc"
#include "isajet/wcon.inc"
C
      REAL AMASS
      REAL AM12
      INTEGER I,J,N,IQ1,IQ2,IW,K
      INTEGER LISTJ(25),LISTW(4)
C
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     $11,-11,12,-12,13,-13,14,-14,15,-15,16,-16/
      DATA LISTW/10,80,-80,90/
C
C          Initialize
C
      IF(IHTYPE.EQ.0) THEN
        WRITE(ITLIS,*) ' YOU MUST SELECT AN HTYPE FOR SUSY HIGGS'
        WRITE(ITLIS,*) ' JOB TERMINATED'
        STOP99
      ENDIF
      HMASS=AMASS(IHTYPE)
      HGAM=0.
      DO 100 I=1,85
        DO 110 J=1,85
          HGAMSS(I,J)=0
110     CONTINUE
100   CONTINUE
C
C          Extract widths from SSMODE common block
C          Note the only 3-body modes are Zff or Wff
C          These are added to the ZZ and WW entries in HCONSS,
C          and the Z* or W* decay is generated later, as for SM Higgs
C
      DO 200 N=1,NSSMOD
        IF(ISSMOD(N).NE.IHTYPE) GO TO 200
        HGAM=HGAM+GSSMOD(N)
        IF(JSSMOD(3,N).NE.0) THEN
C          3-body modes
          IF(IABS(JSSMOD(1,N)).EQ.80) THEN
            HGAMSS(78,79)=HGAMSS(78,79)+0.5*GSSMOD(N)
            HGAMSS(79,78)=HGAMSS(79,78)+0.5*GSSMOD(N)
          ELSEIF(JSSMOD(1,N).EQ.90) THEN
            HGAMSS(80,80)=HGAMSS(80,80)+GSSMOD(N)
          ELSE
            WRITE(ITLIS,1000) ISSMOD(N),(JSSMOD(K,N),K=1,5)
1000        FORMAT(' SETHSS: UNEXPECTED MODE ',I8,' --> ',5I8)
            STOP 99
          ENDIF
          GO TO 200
        ELSE
C          2-body modes
          DO 210 I=1,85
            IF(JSSMOD(1,N).NE.LISTSS(I)) GO TO 210
            DO 220 J=1,85
              IF(JSSMOD(2,N).NE.LISTSS(J)) GO TO 220
              HGAMSS(I,J)=HGAMSS(I,J)+.5*GSSMOD(N)
              HGAMSS(J,I)=HGAMSS(J,I)+.5*GSSMOD(N)
              GO TO 200
220         CONTINUE
210       CONTINUE
        ENDIF
        WRITE(ITLIS,1000) ISSMOD(N),(JSSMOD(K,N),K=1,5)
        STOP99
200   CONTINUE
C
C          W* and Z* mass limits
C
      DO 300 I=1,2
        ZSTARS(1,I)=0.
        DO 310 IW=2,4
          ZSTARS(IW,I)=AMASS(LISTW(IW))
          DO 320 IQ1=2,25
            IQ2=MATCH(IQ1,IW)
            IF(IQ2.EQ.0) GO TO 320
            IF(GOWW(IQ1,I).AND.GOWW(IQ2,I)) THEN
              AM12=AMASS(LISTJ(IQ1))+AMASS(LISTJ(IQ2))+1.0
              ZSTARS(IW,I)=MIN(ZSTARS(IW,I),AM12)
            ENDIF
320       CONTINUE
310     CONTINUE
300   CONTINUE
      RETURN
      END

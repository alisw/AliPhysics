#include "isajet/pilot.h"
      INTEGER FUNCTION ITRANS(ID,MCONV)
C
C...convert (MCONV=1) from ISAJET numbering to PDG numbering
C...     or (MCONV=2) from PDG  numbering to ISAJET numbering
C...
C...Ver 7.21: add extra mesons with IABS(ID) > 10000; these only occur
C...          in a few B decays.
C
C...Thanks to Lynn Garren, Fermilab.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
C
      INTEGER ID,MCONV
      INTEGER I2,I3,J,IDA,IF1,IF2,IND,ITMP,IF3,JS,J1,IS1,IS2,IS3,I4,I1
C
C... ITABI(I) converts miscellaneous ISAJET particle ID's to standard 
C... scheme
      INTEGER ITABI(99,2), NOANT(15)
      SAVE ITABI,NOANT
      DATA ITABI/2,1,3,4,5,6,7,8,21,22,
     1           12,11,14,13,16,15,0,0,0,310,
     2           42,41,43,44,45,46,0,0,47,67,
     3           52,51,54,53,56,55,0,0,77,68,
     4           62,61,63,64,65,66,0,0,78,69,
     5           72,71,74,73,76,75,0,0,0,70,
     6           0,0,0,0,0,0,0,0,0,0,
     7           0,0,0,0,0,0,0,0,0,24,
     8           25,57,58,59,40,33,37,34,38,23,
     9           0,0,0,0,0,0,0,0,0,
     *           2,1,3,4,5,6,7,8,0,0,
     1           12,11,14,13,16,15,0,0,0,0,
     2           9,10,90,80,81,0,0,0,0,0,
     3           0,82,86,88,0,83,87,89,84,85,
     4           22,21,23,24,25,26,29,0,0,0,
     5           32,31,34,33,36,35,82,83,84,0,
     6           42,41,43,44,45,46,30,40,50,60,
     7           52,51,54,53,56,55,39,49,0,0,
     8           0,0,0,0,0,0,0,0,0,0,
     9           0,0,0,0,0,0,0,0,0/
      DATA NOANT/-21,-22,-23,-25,-30,-35,-47,-48,-57,-58,-59,
     1           -67,-68,-69,-70/
C
      IDA=IABS(ID)
      ITRANS=0
      IF(MCONV.NE.1) GO TO 200
C...ISAJET algorithm routine
      CALL FLAVOR(ID,IF1,IF2,IF3,JS,IND)
      IF(IDA.EQ.0) THEN
        WRITE(ITLIS,*) ' ITRANS: particle ID is zero'
      ELSEIF(IDA.LT.100) THEN
        ITRANS=ISIGN(ITABI(IDA,1),ID)
        IF(ID.EQ.-20) ITRANS=130
C...check for illegal antiparticles
        ITMP=ITRANS
        IF(ITMP.LT.0) THEN
          DO 101 J=1,15
            IF(ITMP.EQ.NOANT(J)) ITRANS=0
 101      CONTINUE
        ENDIF
      ELSEIF(IND.NE.0.AND.IDA.LT.10000) THEN
        IS1=IABS(IF1)
        IS2=IABS(IF2)
        IS3=IABS(IF3)
C...mesons
        IF(IS1.EQ.0) THEN
          IF(IS2.LE.2 .AND. IS3.LE.2) THEN
C...         don't change
          ELSE
C...         u and d have opposite definitions
            IF(IS2.LE.2) IS2=ITABI(IS2,1)
            IF(IS3.LE.2) IS3=ITABI(IS3,1)
          ENDIF
          ITRANS=IS3*100 + IS2*10 + 2*JS+1
          ITRANS=ISIGN(ITRANS,ID)
C...         charmed and top mesons have wrong sign
          IF(IS3.EQ.4 .AND. IS2.NE.4) ITRANS=-ITRANS
          IF(IS3.EQ.6 .AND. IS2.NE.6 .AND. IS2.NE.4) ITRANS=-ITRANS
C...check for illegal antiparticles
          IF(IS2.EQ.IS3 .AND. ID.LT.0) ITRANS=0
C...diquarks
        ELSEIF(IS3.EQ.0) THEN
C...         u and d have opposite definitions
          IF(IS1.LE.2) IS1=ITABI(IS1,1)
          IF(IS2.LE.2) IS2=ITABI(IS2,1)
          IF(IS2.LT.IS1) THEN
            ITRANS=IS1*1000 + IS2*100 + 2*JS+1
          ELSE
            ITRANS=IS2*1000 + IS1*100 + 2*JS+1
          ENDIF
          ITRANS=ISIGN(ITRANS,ID)
C...         charmed and top mesons have wrong sign
          IF(IS2.EQ.4 .AND. IS1.NE.4) ITRANS=-ITRANS
          IF(IS2.EQ.6 .AND. IS1.NE.6 .AND. IS1.NE.4) ITRANS=-ITRANS
C...baryons
        ELSE
C...         u and d have opposite definitions
          IF(IS1.LE.2) IS1=ITABI(IS1,1)
          IF(IS2.LE.2) IS2=ITABI(IS2,1)
          IF(IS3.LE.2) IS3=ITABI(IS3,1)
          IF(IS3.LE.2) THEN
            ITRANS=IS1*1000 + IS2*100 + IS3*10 + 2*JS+2
          ELSEIF(IS1.LE.2 .AND. IS2.LE.2) THEN
            ITRANS=IS3*1000 + IS1*100 + IS2*10 + 2*JS+2
          ELSE
            ITRANS=IS3*1000 + IS2*100 + IS1*10 + 2*JS+2
          ENDIF
          ITRANS=ISIGN(ITRANS,ID)
        ENDIF
      ELSEIF(IND.GT.0.AND.IDA.GT.10000) THEN
C...Special mesons. 
        IF(ID.EQ.10121) THEN
          ITRANS=20213
        ELSEIF(ID.EQ.10111) THEN
          ITRANS=20113
        ELSEIF(ID.EQ.10131) THEN
          ITRANS=10323
        ELSEIF(ID.EQ.10231) THEN
          ITRANS=10313
        ELSEIF(ID.EQ.30131) THEN
          ITRANS=30323
        ELSEIF(ID.EQ.30231) THEN
          ITRANS=30313
        ELSEIF(ID.EQ.10110) THEN
          ITRANS=10221
        ELSEIF(ID.EQ.10441) THEN
          ITRANS=20443
        ELSEIF(ID.EQ.20440) THEN
          ITRANS=10441
        ELSEIF(ID.EQ.20441) THEN
          ITRANS=10443
        ELSEIF(ID.EQ.20442) THEN
          ITRANS=445
        ELSE
          ITRANS=ID
        ENDIF
      ENDIF
      GO TO 300
C
 200  IF(MCONV.NE.2) GO TO 300
      J1=MOD(IDA,10)
      I1=MOD(IDA/10,10)
      I2=MOD(IDA/100,10)
      I3=MOD(IDA/1000,10)
      I4=MOD(IDA/10000,10)
      IF(IDA.EQ.0) THEN
        WRITE(ITLIS,*) ' ITRANS: particle ID is zero'
C...elementary particles
      ELSEIF(IDA.LT.100) THEN
        ITRANS=ISIGN(ITABI(IDA,2),ID)
C...check for illegal antiparticles
        IF(ID.LT.0) THEN
          DO 201 J=1,15
            IF(ID.EQ.NOANT(J)) ITRANS=0
 201      CONTINUE
        ENDIF
C...K short and K long
      ELSEIF(ID.EQ.130) THEN
        ITRANS=-20
      ELSEIF(ID.EQ.310) THEN
        ITRANS=20
C...mesons
      ELSEIF(I3.EQ.0) THEN
          IF(I1.LE.2 .AND. I2.LE.2) THEN
C...         don't change
          ELSE
C...         u and d have opposite definitions
            IF(I1.LE.2) I1=ITABI(I1,2)
            IF(I2.LE.2) I2=ITABI(I2,2)
          ENDIF
          ITRANS=I1*100 + I2*10 + (J1-1)/2
          ITRANS=ISIGN(ITRANS,ID)
C...         charmed and top mesons have wrong sign
          IF(I2.EQ.4 .AND. I1.NE.4) ITRANS=-ITRANS
          IF(I2.EQ.6 .AND. I1.NE.6 .AND. I1.NE.4) ITRANS=-ITRANS
C...check for illegal antiparticles
          IF(I2.EQ.I1 .AND. ID.LT.0) ITRANS=0
C...diquarks
      ELSEIF(I1.EQ.0) THEN
C...         u and d have opposite definitions
          IF(I3.LE.2) I3=ITABI(I3,2)
          IF(I2.LE.2) I2=ITABI(I2,2)
          IF(I3.LT.I2) THEN
            ITRANS=I3*1000 + I2*100 + (J1-1)/2
          ELSE
            ITRANS=I2*1000 + I3*100 + (J1-1)/2
          ENDIF
          ITRANS=ISIGN(ITRANS,ID)
C...         charmed and top mesons have wrong sign
          IF(I2.EQ.4 .AND. I3.NE.4) ITRANS=-ITRANS
          IF(I2.EQ.6 .AND. I3.NE.6 .AND. I3.NE.4) ITRANS=-ITRANS
C...baryons
      ELSE
C...         u and d have opposite definitions
          IF(I3.LE.2) I3=ITABI(I3,2)
          IF(I2.LE.2) I2=ITABI(I2,2)
          IF(I1.LE.2) I1=ITABI(I1,2)
          IF(I3.LE.2) THEN
            ITRANS=I3*1000 + I2*100 + I1*10 + (J1-2)/2
          ELSEIF(I1.LE.2 .AND. I2.LE.2) THEN
            ITRANS=I2*1000 + I1*100 + I3*10 + (J1-2)/2
          ELSE
            ITRANS=I1*1000 + I2*100 + I3*10 + (J1-2)/2
          ENDIF
          ITRANS=ISIGN(ITRANS,ID)
      ENDIF

 300  RETURN
      END

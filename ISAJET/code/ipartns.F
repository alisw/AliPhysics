#include "isajet/pilot.h"
      SUBROUTINE IPARTNS(NPRTNS,IDS,PRTNS,IDQ,WEIGHT,WZDK)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     fill PJETS array from a list of input partons
C-   Inputs  : 
C-     NPRTNS          = number of partons
C-     IDS(NPRTNS)     = parton ids
C-     PRTNS(4,NPRTNS) = parton 4 vectors
C-     IDQ(2)          = initial partons
C-     WEIGHT          = weight
C-     WZDK            = if true last 2 partons are from W,Z decay
C-     
C-
C-   Created   8-OCT-1991   Serban D. Protopopescu
C-   Updated  17-APR-1996   Serban D. Protopopescu  
C-    added entry evcuts to supply evolution limits
C-    modified DrellYan (keys(3)) to stay within VECBOS jet ranking 
C-   Updated  16-JUN-1998   F. Paige
C-    Removed ISAZEB dependence: use ISPJET and do not call ISPETA
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      INTEGER NPRTNS,IDS(NPRTNS),IDQ(2)
      REAL    PRTNS(4,NPRTNS),WEIGHT
      LOGICAL WZDK
#include "isajet/final.inc"
#include "isajet/idrun.inc"
#include "isajet/jetpar.inc"
#include "isajet/keys.inc"
#include "isajet/nodcay.inc"
#include "isajet/partcl.inc"
#include "isajet/pjets.inc"
#include "isajet/primar.inc"
#include "isajet/q1q2.inc"
#include "isajet/totals.inc"
      REAL    SUM(4),AMASS
      INTEGER K,J,IWZ,ID,NQS
      INTEGER MAXQ
      PARAMETER (MAXQ=15)
      INTEGER I,NP,JDORD(MAXQ),JIORD(MAXQ),NPJ
      REAL    ETAQ(MAXQ),PHIQ(MAXQ),THQ(MAXQ),PTQ(MAXQ)
      REAL    ETCUT,ETIN,RCUT,RIN,R
      REAL    PPI
      REAL    PXPT(MAXQ),PXETA(MAXQ),PXPHI(MAXQ)
      LOGICAL DOEVOL,DOEVIN
      DOUBLE PRECISION PI, TWOPI, HALFPI, RADIAN
      PARAMETER (PI=        3.1415 92653 58979 32384 6 D0)
      PARAMETER (TWOPI=     6.2831 85307 17958 64769 3 D0)
      PARAMETER (HALFPI=    1.5707 96326 79489 66192 3 D0)
      PARAMETER (RADIAN= 0.0174532 92519 94329 5769237 D0)
C----------------------------------------------------------------------
C
      NJET=0
C
C          handle W's and Z's
C          
      IEVT=IEVT+1 
      IWZ=0
      NQS=NPRTNS
      IF(WZDK) NQS=NPRTNS-2
      DO 1 J=1,NPRTNS
        ID=IABS(IDS(J))
        IF(ID.GT.79) THEN
          IF(ID.EQ.90) JWTYP=4
          IF(IDS(J).EQ.80) JWTYP=2
          IF(IDS(J).EQ.-80) JWTYP=3
          IDENTW=IDS(J)
          DO 2 K=1,4
            QWJET(K)=PRTNS(K,J)
   2      CONTINUE
          QWJET(5)=SQRT(QWJET(4)**2-QWJET(1)**2-QWJET(2)**2-QWJET(3)**2)
          IWZ=J
        ENDIF
   1  CONTINUE
      DO 4 J=NQS+1,NPRTNS
        ID=IABS(IDS(J))
        NJET=NJET+1
        DO 3 K=1,4
          PJETS(K,NJET)=PRTNS(K,J)
    3   CONTINUE
        IDJETS(NJET)=IDS(J)
        PJETS(5,NJET)=AMASS(ID)
    4 CONTINUE
C          W,Z decays were not in input
      IF(IWZ.NE.0.AND.NJET.EQ.0) THEN    
        NJET=2
        CALL ISWDKY
      ENDIF
C
C      fill with the other partons
C
      DO 5 K=1,4
        SUM(K)=0
   5  CONTINUE
      DO 11 J=1,NQS
        ID=IABS(IDS(J))
        IF(IWZ.NE.J.AND.ID.LT.11) THEN
          NJET=NJET+1
          IDJETS(NJET)=IDS(J)
          DO 12 K=1,4
            PJETS(K,NJET)=PRTNS(K,J)
  12      CONTINUE
          PJETS(5,NJET)=PRTNS(4,J)**2-PRTNS(1,J)**2-PRTNS(2,J)**2-
     $      PRTNS(3,J)**2
          IF ( PJETS(5,NJET).GT.0. ) THEN
            PJETS(5,NJET)=SQRT(PJETS(5,NJET))
          ELSE
            PJETS(4,NJET)=SQRT(PRTNS(4,J)**2-PJETS(5,NJET))
            PJETS(5,NJET)=0.
          ENDIF
        ENDIF
        DO 13 K=1,4
          SUM(K)=SUM(K)+PRTNS(K,J)
  13    CONTINUE
  11  CONTINUE
C
C        eta and phi of incoming partons 
      IF(DOEVOL) THEN
        NP=NQS-1
        DO 114 I=1,NP
          PPI=SQRT(PRTNS(1,I)**2+PRTNS(2,I)**2+PRTNS(3,I)**2)
          IF(PPI.GT.0.AND.PPI.GT.ABS(PRTNS(3,I))) THEN
            THQ(I)=ACOS(PRTNS(3,I)/PPI)
            ETAQ(I)=-LOG(TAN(THQ(I)/2))
          ELSE
            THQ(I)=0
            ETAQ(I)=SIGN(999.,PRTNS(3,I))
          ENDIF
          PTQ(I)=SQRT(PRTNS(1,I)**2+PRTNS(2,I)**2)
          IF(PTQ(I).GT.0) THEN
            PHIQ(I)=ATAN2(PRTNS(2,I),PRTNS(1,I))
            IF(PHIQ(I).LT.0) PHIQ(I)=PHIQ(I)+TWOPI
          ELSE
            PHIQ(I)=0
          ENDIF
 114    CONTINUE
C
C ... Order partons in pt
C
        DO 115 I = 1 , NP
          JIORD(I) = I
          PXPT(I)=PTQ(I)
 115    CONTINUE
        CALL ISASRT(PXPT(1),NP,JIORD)
        DO 116 I = 1 , NP
          PXPT(I)=PTQ(I)
          PXETA(I)=ETAQ(I)
          PXPHI(I)=PHIQ(I)
          JDORD(I) = JIORD(NP-I+1)
 116    CONTINUE
        DO 117 I = 1 , NP
          PTQ(I)=PXPT(JDORD(I))
          ETAQ(I)=PXETA(JDORD(I))
          PHIQ(I)=PXPHI(JDORD(I))
 117    CONTINUE
      ENDIF
C
C
  15  CONTINUE
      PBEAM(1)=(ECM-SUM(4)-SUM(3))/2.
      PBEAM(2)=(ECM-SUM(4)+SUM(3))/2.
      QSQ=SQRT(SUM(4)**2-SUM(3)**2-SUM(2)**2-SUM(1)**2)
      CALL RANFMT 
      NPTCL=0
      IF(KEYS(3)) THEN
        STDDY=.FALSE.
        IF(NQS.EQ.1.OR.NJET.LT.3) STDDY=.TRUE.
      ENDIF
      CALL IPRTNS(NQS,PRTNS,IDQ)
      IF(.NOT.NOEVOL) THEN
        CALL EVOLVE
C
C            special check for VECBOS
        IF(DOEVOL) THEN  
C       Find parton jets
          CALL ISPJET(RCUT,ETCUT,NPJ,PXPT,PXPHI,PXETA)  
          IF(NPJ.GE.NP.AND.PXPT(NP).GT.PTQ(NP)) THEN
            R=SQRT((PXETA(NP)-ETAQ(NP))**2+(PXPHI(NP)-PHIQ(NP))**2)
            IF(R.GT.RCUT) GOTO 15
          ENDIF
        ENDIF
C
        IF(.NOT.NOHADR) THEN
          CALL FRGMNT
          CALL MBIAS
        ENDIF
      ENDIF
      WT=WEIGHT
      SUMWT=SUMWT+WT
      SIGF=SUMWT
      NKINF=IEVT
      NEVENT=IEVT
  999 RETURN
C
C     Entry point to set parameters
C
      ENTRY EVCUTS(RIN,ETIN,DOEVIN)
      RCUT=RIN
      ETCUT=ETIN
      DOEVOL=DOEVIN
      RETURN
      END

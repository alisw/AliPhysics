#include "isajet/pilot.h"
      SUBROUTINE ISPJET(DRCUT,ETCUT,NPJ,PJPT,PJPHI,PJETA)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : COMBINES PARTONS INTO PARTON JETS
C-                         based on PJCONE
C-   Inputs
C-        DRCUT  - dR=sqrt(dETA**2+dPHI**2) cut around Leading Partons.
C-        ETCUT  - Transverse Energy cut (minimum for defining a JET ).
C-
C-   Outputs : 
C-        NPJ       = No. of Parton Jets found.
C-        PJPT(NPJ) = pt of partons
C-        PJPHI(NPJ)= phi "
C-        PJETA(NPJ)= eta "
C-        
C-   created  16-APR-1996   Serban D. Protopopescu   
C-   Updated  16-JUN-1998   F. Paige
C-    Copy of ISAZEB routine ISA_PJETS to be used by IPARTNS
C-    
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/pjets.inc"
#include "isajet/jetset.inc"
C
      INTEGER NPJ
      REAL    DRCUT
      REAL    PJPHI(*), PJETA(*), PJPT(*)
      INTEGER NP, JP, J, JO, JOP1, JOP2, JP1, JP2, ISKP, IP
      REAL    X1, Y1, PHI1, PHI2,TH
      REAL    DETA, DPHI, DR,  ETCUT
      INTEGER NPMAX
      PARAMETER (NPMAX=50)
      INTEGER JIORD(NPMAX), JDORD(NPMAX), JCNN(NPMAX,NPMAX)
      INTEGER JSKP(NPMAX)
      INTEGER I, JJ, K
      REAL    PJIN(4,NPMAX), PINPHI(NPMAX), PINETA(NPMAX)
      REAL    PINPT(NPMAX),PDMPT(NPMAX)
      REAL    PJ(4,NPMAX)
      REAL    EPS
      DOUBLE PRECISION PI, TWOPI, HALFPI, RADIAN
C
C last significant (correctly rounded) decimal place on VAX:
C                                                |
C                                                V
        PARAMETER (PI=        3.1415 92653 58979 32384 6 D0)
        PARAMETER (TWOPI=     6.2831 85307 17958 64769 3 D0)
        PARAMETER (HALFPI=    1.5707 96326 79489 66192 3 D0)
        PARAMETER (RADIAN= 0.0174532 92519 94329 5769237 D0)
C
      PARAMETER( EPS = 1.0E-5 )
C----------------------------------------------------------------------
C
      NP=0
      DO 10 I=1,NJSET
        IF(JDCAY(I).EQ.0.AND.IABS(JTYPE(I)).LT.10) THEN
          NP = NP + 1
          DO 11 K=1,4
            PJIN(K,NP)=PJSET(K,I)
   11     CONTINUE
          PINPT(NP)  = SQRT( PJIN(1,NP)**2+PJIN(2,NP)**2 )
          PINPHI(NP) = ATAN2 (PJIN(2,NP),PJIN(1,NP)+EPS)
          IF(PINPHI(NP).LT.0.)PINPHI(NP)=PINPHI(NP)+TWOPI
          TH = ATAN2 (PINPT(NP),PJIN(3,NP)+EPS)
          PINETA(NP) = -ALOG ( ABS(TAN(TH/2.)) + EPS )
          IF(NP.GE.NPMAX) GOTO 35
        ENDIF
   10 CONTINUE
   35 CONTINUE    ! jump here if more than NPMAX partons
C
C ... Order partons in pt
C
      DO 100 JP = 1 , NP
        JIORD(JP) = JP
  100 PDMPT(JP)=PINPT(JP)
      CALL ISASRT(PDMPT(1),NP,JIORD)
      DO 105 JP = 1 , NP
  105 JDORD(JP) = JIORD(NP-JP+1)
C
C ... Combine partons close in r space
C
      DO 110 J = 1 , NP
        JO=JDORD(J)
  110 JCNN(JO,1)=0
      ISKP=0
      DO 120 JP1 = 1 , NP-1
        JOP1=JDORD(JP1)
C ... Check if parton already connected to other one
        IF ( JCNN(JOP1,1).EQ.-1 ) GOTO 120
        DO 130 JP2 = JP1+1 , NP
          JOP2=JDORD(JP2)
C ... Check if parton already connected to other one
          IF ( JCNN(JOP2,1).EQ.-1 ) GOTO 130
          DETA = PINETA(JOP1) - PINETA(JOP2)
          PHI1 = PINPHI(JOP1)
          PHI2 = PINPHI(JOP2)
          X1 = COS(PHI2-PHI1)
          Y1 = SIN(PHI2-PHI1)
          IF(X1.EQ.0.0) THEN
            DPHI = HALFPI
          ELSE
            DPHI = ATAN2(Y1,X1)
          END IF
          DR = SQRT(DETA**2+DPHI**2)
C --- Criterion for combining partons
          IF ( DR.LT.DRCUT ) THEN
            JCNN(JOP1,1)=JCNN(JOP1,1)+1
            JCNN(JOP2,1)=-1
            JCNN(JOP1,JCNN(JOP1,1)+1)=JOP2
            ISKP=ISKP+JCNN(JOP1,1)
            JSKP(ISKP)=JOP2
          ELSE
            GOTO 130
          ENDIF
  130   CONTINUE
  120 CONTINUE
C
C ... Bookkeeping for parton jets
C
      DO 150 IP = 1 , NPJ
        PJPHI(IP)=0.
        PJETA(IP)=0.
        PJPT(IP) =0.
  150 CONTINUE
      NPJ=0
      DO 200 JP1 = 1 , NP
        JOP1=JDORD(JP1)
C ... Already connected, single parton, or has others to connect to
        IF ( JCNN(JOP1,1).GE.0 ) THEN
          NPJ=NPJ+1
          DO 151 K=1,4
            PJ(K,NPJ)=PJIN(K,JOP1)
  151     CONTINUE
          PJPHI(NPJ) = PINPHI(JOP1)
          PJETA(NPJ) = PINETA(JOP1)
          PJPT(NPJ) = PINPT(JOP1)
          IF ( JCNN(JOP1,1).EQ.0 ) GOTO 205
          DO 210 JJ = 1 , JCNN(JOP1,1)
            PJ(1,NPJ) = PJ(1,NPJ) + PJIN(1,JCNN(JOP1,JJ+1))
            PJ(2,NPJ) = PJ(2,NPJ) + PJIN(2,JCNN(JOP1,JJ+1))
            PJ(3,NPJ) = PJ(3,NPJ) + PJIN(3,JCNN(JOP1,JJ+1))
            PJ(4,NPJ) = PJ(4,NPJ) + PJIN(4,JCNN(JOP1,JJ+1))
  210     CONTINUE
          PJPT(NPJ) = SQRT( PJ(1,NPJ)**2 + PJ(2,NPJ)**2 )
          PJPHI(NPJ) = ATAN2 (PJ(2,NPJ),PJ(1,NPJ)+EPS)
          IF(PJPHI(NPJ).LT.0.)PJPHI(NPJ)=PJPHI(NPJ)+TWOPI
          TH = ATAN2 (PJPT(NPJ),PJ(3,NPJ)+EPS)
          PJETA(NPJ) = -ALOG ( ABS(TAN(TH/2.)) + EPS )
C ... Criterion for dropping a parton jet ( et < etcut )
  205     IF ( PJPT(NPJ).GT.ETCUT ) GOTO 200
          NPJ=NPJ-1
        ENDIF
  200 CONTINUE
C
C ... Order pjets in pt
C
      DO 300 JP = 1 , NPJ
        JIORD(JP) = JP
  300 PDMPT(JP)=PJPT(JP)
      CALL ISASRT(PDMPT(1),NPJ,JIORD)
      DO 305 JP = 1 , NPJ
        PINPT(JP)=PJPT(JP)
        PINETA(JP)=PJETA(JP)
        PINPHI(JP)=PJPHI(JP)
  305 JDORD(JP) = JIORD(NPJ-JP+1)
      DO 306 JP = 1 , NPJ
        PJPT(JP)=PINPT(JDORD(JP))
        PJETA(JP)=PINETA(JDORD(JP))
        PJPHI(JP)=PINPHI(JDORD(JP))
  306 CONTINUE
C-
  999 RETURN
      END

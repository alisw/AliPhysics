#include "isajet/pilot.h"
      SUBROUTINE MULJET(WT)
C
C          Using masses from /MGKIN/, generate NJET<=MXJETS body phase 
C          space point satisfying cuts:
C          (1) Generate kinematic point using successive 2-body decays
C              with Jacobean
C              dPhi_N(p1...pN) = dQ1 p1/(4*pi) dPhi_(N-1)(q1...pN)
C          (2) Apply individual jet cuts from /JETLIM/ and dijet
C              cuts from /MGLIMS/ to ensure IR-safe cross section.
C          (3) Return weight WT or 0 if outside limits. 
C
C          Note that WT contains various constant factors that were 
C          dropped in DECAY:
C            1/(2*SHMG) Jacobean
C            Jacobean for dQ = (EHMG-SUM)*dRANF
C            Factors of 4pi
C
C          MadGraph/Helas notation: 
C          PJETS8(0:3,1:2) = initial momenta
C          PJETS8(0:3,3:NJET+2) = final momenta
C          Note: ANSI extensions, e.g. REAL*8 P(0:3) are required for 
C          compatibility with Helas and MadGraph. :-(
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
C
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/mglims.inc"
#include "isajet/pjets.inc"
#include "isajet/mgkin.inc"
#include "isajet/primar.inc"
C
C          Local variables; MXJETS defined in /PJETS/
C
      REAL*8 PGEN(0:3,MXJETS),AMGEN(MXJETS),RND(MXJETS)
      REAL*8 SHMG,EHMG,YHMG,SUM,SUM1,RNEW,WT,QCM,PI,
     $U(3),PHI,BETA(3),GAMMA,BP,PTI,PPI,YI,XJI,PHII,AMIJ
      REAL*8 CYHMG,SYHMG,E1,E2,P12,DELTAQ
      REAL*8 PCM,A,B,C
      INTEGER NJET1,I,JJ1,NTRY,J,JSAVE,II,K
      REAL RANF
C
C          Function definition
C
      PCM(A,B,C)=SQRT((A-B-C)*(A+B+C)*(A-B+C)*(A+B-C))/(2.*A)
C
C
C          Generate COM mass and rapidity
C
      PI=4.D0*DATAN(1.D0)
100   CONTINUE
      SHMG=EHMGMN**2+(EHMGMX**2-EHMGMN**2)*RANF()
      EHMG=SQRT(SHMG)
      YHMG=YHMGMN+(YHMGMX-YHMGMN)*RANF()
      IF(EHMG*EXP(ABS(YHMG)).GT.ECM) GO TO 999
      IF(EHMG.LT.AMJET8(1)+AMJET8(2)) GO TO 999
      CYHMG=DCOSH(YHMG)
      SYHMG=SINH(YHMG)
      AMGEN(1)=EHMG
      PGEN(1,1)=0
      PGEN(2,1)=0
      PGEN(3,1)=EHMG*SYHMG
      PGEN(0,1)=EHMG*CYHMG
      E1=(EHMG**2+AMJET8(1)**2-AMJET8(2)**2)/(2*EHMG)
      E2=(EHMG**2-AMJET8(1)**2+AMJET8(2)**2)/(2*EHMG)
      P12=PCM(EHMG,AMJET8(1),AMJET8(2))
C          Initial momenta
      PJETS8(1,1)=0
      PJETS8(2,1)=0
      PJETS8(3,1)=SYHMG*E1+CYHMG*P12
      PJETS8(0,1)=CYHMG*E1+SYHMG*P12
      PJETS8(1,2)=0
      PJETS8(2,2)=0
      PJETS8(3,2)=SYHMG*E1-CYHMG*P12
      PJETS8(0,2)=CYHMG*E1-SYHMG*P12
C
      NJET1=NJET-1
      SUM=0
      DO 110 I=1,NJET
        SUM=SUM+AMJET8(I+2)
110   CONTINUE
      IF(SUM.GE.EHMG) GO TO 999
      DELTAQ=EHMG-SUM
C
C          Generate masses for uniform NJET-body phase space.
C
      NTRY=0
200   CONTINUE
      NTRY=NTRY+1
      IF(NTRY.GT.NTRIES) THEN
        WRITE(ITLIS,9999) NTRY
9999    FORMAT(//2X,'ERROR IN MULJET ... NTRY = ',I8)
        STOP99
      ENDIF
      RND(1)=1
      DO 210 I=2,NJET1
        RNEW=RANF()
        DO 220 JJ1=1,I-1
          J=I-JJ1
          JSAVE=J+1
          IF(RNEW.LE.RND(J)) GO TO 210
          RND(JSAVE)=RND(J)
220     CONTINUE
210   RND(JSAVE)=RNEW
      RND(NJET)=0
C          Jacobean for d(shmg)d(yhmg) and overall 1/(2*shmg)
      WT=(EHMGMX**2-EHMGMN**2)*(YHMGMX-YHMGMN)/(2*SHMG)
      SUM1=SUM
      DO 230 I=2,NJET
        SUM1=SUM1-AMJET8(I-1+2)
        AMGEN(I)=SUM1+RND(I)*(AMGEN(1)-SUM)
        IF(AMGEN(I-1).LE.AMGEN(I)+AMJET8(I-1+2)) GO TO 200
C          Jacobean for sigma_n -> sigma_n-1
        WT=WT*PCM(AMGEN(I-1),AMGEN(I),AMJET8(I-1+2))*DELTAQ/(4*PI**2)
230   CONTINUE
C          Jacobean for final 2-body decay differs by this factor
      WT=WT*PI/(DELTAQ*EHMG)
C
C          Carry out 2-body decays
C
      DO 310 I=1,NJET1
        QCM=PCM(AMGEN(I),AMGEN(I+1),AMJET8(I+2))
        U(3)=2.*RANF()-1
        PHI=2*PI*RANF()
        U(1)=SQRT(1-U(3)**2)*COS(PHI)
        U(2)=SQRT(1-U(3)**2)*SIN(PHI)
        DO 320 J=1,3
          PJETS8(J,I+2)=QCM*U(J)
          PGEN(J,I+1)=-PJETS8(J,I+2)
320     CONTINUE
        PJETS8(0,I+2)=SQRT(QCM**2+AMJET8(I+2)**2)
        PGEN(0,I+1)=SQRT(QCM**2+AMGEN(I+1)**2)
310   CONTINUE
C
      DO 330 J=0,3
        PJETS8(J,NJET+2)=PGEN(J,NJET)
330   CONTINUE
C
C          Boost PGEN frames to lab frame.
C
      DO 400 II=1,NJET1
        I=NJET-II
        DO 410 J=1,3
          BETA(J)=PGEN(J,I)/PGEN(0,I)
410     CONTINUE
        GAMMA=PGEN(0,I)/AMGEN(I)
        DO 420 K=I,NJET
          BP=BETA(1)*PJETS8(1,K+2)+BETA(2)*PJETS8(2,K+2)+
     $    BETA(3)*PJETS8(3,K+2)
          DO 430 J=1,3
            PJETS8(J,K+2)=PJETS8(J,K+2)+GAMMA*BETA(J)*(PJETS8(0,K+2)
     $      +BP*GAMMA/(GAMMA+1.))
430       CONTINUE
          PJETS8(0,K+2)=GAMMA*(PJETS8(0,K+2)+BP)
420     CONTINUE
400   CONTINUE
C
C          Check limits
C
      DO 500 I=1,NJET
        PTI=SQRT(PJETS8(1,I+2)**2+PJETS8(2,I+2)**2)
        IF(PTI.LE.PTMIN(I).OR.PTI.GE.PTMAX(I)) GO TO 999
        PPI=SQRT(PTI**2+PJETS8(3,I+2)**2)
        IF(PPI.LE.PMIN(I).OR.PPI.GE.PMAX(I)) GO TO 999
        XJI=PJETS8(3,I+2)/PPI
        IF(XJI.LE.XJMIN(I).OR.XJI.GE.XJMAX(I)) GO TO 999
        PHII=ATAN2(PJETS8(2,I+2),PJETS8(1,I+2))
        IF(PHII.LT.0) PHII=PHII+2*PI
        IF(PHII.LE.PHIMIN(I).OR.PHII.GE.PHIMAX(I)) GO TO 999
        YI=-LOG(TAN(ACOS(XJI)/2))
        IF(YI.LE.YJMIN(I).OR.YI.GE.YJMAX(I)) GO TO 999
500   CONTINUE
C
      DO 510 I=1,NJET
        DO 520 J=I+1,NJET
          AMIJ=(PJETS8(0,I+2)+PJETS8(0,J+2))**2
     $        -(PJETS8(1,I+2)+PJETS8(1,J+2))**2
     $        -(PJETS8(2,I+2)+PJETS8(2,J+2))**2
     $        -(PJETS8(3,I+2)+PJETS8(3,J+2))**2
          AMIJ=SIGN(SQRT(ABS(AMIJ)),AMIJ)
          IF(AMIJ.LE.AMIJMN(I,J).OR.AMIJ.GE.AMIJMX(I,J)) GO TO 999
520     CONTINUE
510   CONTINUE
C
      RETURN
C
999   WT=0
      RETURN
      END

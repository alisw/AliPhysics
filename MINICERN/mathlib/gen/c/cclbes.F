*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:01:56  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_SINGLE)
      SUBROUTINE CCLBES(ZZ,ETA1,ZLMIN,NL,FC,GC,FCP,GCP,SIG,KFN,MODE1,
     1                  IFAIL,IPR)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C  Original title : COULCC                                             C
C                                                                      C
C  A. R. Barnett           Manchester  March   1981                    C
C  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions C
C                                                                      C
C  The FORM (not the SUBSTANCE) of this program has been modified      C
C   by K.S. KOLBIG (CERN)    December 1987                             C
C                                                                      C
C  original program  RCWFN       in    CPC  8 (1974) 377-395           C
C                 +  RCWFF       in    CPC 11 (1976) 141-142           C
C                 +  COULFG      in    CPC 27 (1982) 147-166           C
C  description of real algorithm in    CPC 21 (1981) 297-314           C
C  description of complex algorithm    JCP 64 (1986) 490-509           C
C  this version written up       in    CPC 36 (1985) 363-372           C
C                                                                      C
C  CCLBES returns F,G,G',G',SIG for complex ETA, ZZ, and ZLMIN,        C
C   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL inclusive,   C
C   thus giving  complex-energy solutions to the Coulomb Schrodinger   C
C   equation,to the Klein-Gordon equation and to suitable forms of     C
C   the Dirac equation ,also spherical & cylindrical Bessel equations  C
C                                                                      C
C  if ABS(MODE1)                                                       C
C            = 1  get F,G,F',G'   for integer-spaced lambda values     C
C            = 2      F,G      unused arrays must be dimensioned in    C
C            = 3      F,  F'          call to at least length (0:1)    C
C            = 4      F                                                C
C            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )    C
C            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in C
C            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC C
C            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )    C
C                                                                      C
C     if MODE1 < 0 then the values returned are scaled by an exponen-  C
C                  tial factor (dependent only on ZZ) to bring nearer  C
C                  unity the functions for large ABS(ZZ), small ETA ,  C
C                  and ABS(ZL) < ABS(ZZ).                              C
C        Define SCALE = (  0        if MODE1 > 0                       C
C                       (  IMAG(ZZ) if MODE1 < 0  &  KFN < 3           C
C                       (  REAL(ZZ) if MODE1 < 0  &  KFN = 3           C
C        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)                 C
C         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )                   C
C               or EXP(SCALE)       * ( H+, H(1), or K)                C
C               or EXP(-SCALE)      * ( H- or H(2) )                   C
C                                                                      C
C  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G    C
C           =  1   spherical Bessel      "      "     "       j & y    C
C           =  2 cylindrical Bessel      "      "     "       J & Y    C
C           =  3 modified cyl. Bessel    "      "     "       I & K    C
C                                                                      C
C          and where Coulomb phase shifts put in SIG if KFN=0 (not -1) C
C                                                                      C
C  The use of MODE and KFN is independent                              C
C    (except that for KFN=3,  H(1) & H(2) are not given)               C
C                                                                      C
C  With negative orders lambda, CCLBES can still be used but with      C
C    reduced accuracy as CF1 becomes unstable. The user is thus        C
C    strongly advised to use reflection formulae based on              C
C    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi) C
C                                                                      C
C  Precision:  results to within 2-3 decimals of 'machine accuracy',   C
C              except in the following cases:                          C
C              (1) if CF1A fails because X too small or ETA too large  C
C               the F solution  is less accurate if it decreases with  C
C               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)   C
C              (2) if ETA is large (e.g. >> 50) and X inside the       C
C                turning point, then progressively less accuracy       C
C              (3) if ZLMIN is around sqrt(ACCUR) distance from an     C
C               integral order and abs(X) << 1, then errors present.   C
C               RERR traces the main roundoff errors.                  C
C                                                                      C
C   CCLBES is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14  C
C          with a section of doubled REAL*16 for less roundoff errors. C
C          (If no doubled precision available, increase JMAX to eg 100)C
C   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32  C
C   For single precision CDC (48 bits) reassign                        C
C        DOUBLE PRECISION=REAL etc.                                    C
C                                                                      C
C   IPR    on input   = 0 : no printing of error messages              C
C                    ne 0 : print error messages on file 6             C
C   IFAIL  in output = -2 : argument out of range                      C
C                    = -1 : one of the continued fractions failed,     C
C                           or arithmetic check before final recursion C
C                    =  0 : All Calculations satisfactory              C
C                    ge 0 : results available for orders up to & at    C
C                             position NL-IFAIL in the output arrays.  C
C                    = -3 : values at ZLMIN not found as over/underflowC
C                    = -4 : roundoff errors make results meaningless   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     Machine dependent constants :                                    C
C                                                                      C
C     ACCU     target bound on relative error (except near 0 crossings)C
C               (ACCUR should be at least 100 * ACC8)                  C
C     ACC8    smallest number with 1+ACC8 .ne.1 in REAL*8  arithmetic  C
C     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic C
C     FPMAX    magnitude of largest floating point number * ACC8       C
C     FPMIN    magnitude of smallest floating point number / ACC8      C
C                                                                      C
C     Parameters determining region of calculations :                  C
C                                                                      C
C        R20      estimate of (2F0 iterations)/(CF2 iterations)        C
C        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily     C
C        XNEAR    minimum ABS(X) for CF2 to converge accurately        C
C        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series  C
C        JMAX     size of work arrays for Pade accelerations           C
C        NDROP    number of successive decrements to define instabilityC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     C309R1 = CF2,   C309R2 = F11,    C309R3 = F20,
C     C309R4 = CF1R,  C309R5 = CF1C,   C309R6 = CF1A,
C     C309R7 = RCF,   C309R8 = CTIDY.
C
      IMPLICIT COMPLEX (A-H,O-Z)
C
      COMMON /COULC2/ NFP,N11,NPQ1,NPQ2,N20,KAS1,KAS2
      EQUIVALENCE (NPQ(1),NPQ1),(NPQ(2),NPQ2)
      EQUIVALENCE (KAS(1),KAS1),(KAS(2),KAS2)
      REAL ZERO,ONE,TWO,HALF
      REAL R20,ASYM,XNEAR
      REAL FPMAX,FPMIN,FPLMAX,FPLMIN
      REAL ACCU,ACC8,ACC16
      REAL HPI,TLOG
      REAL ERR,RERR,ABSC,ACCUR,ACCT,ACCH,ACCB,C309R4
      REAL PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,ABSX
      REAL SX1,SETA,SZLL

      LOGICAL LPR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2,NPINT

      PARAMETER(ZERO = 0, ONE = 1, TWO = 2, HALF = ONE/TWO)
      PARAMETER(CI = (0,1), CIH = CI/TWO)
      PARAMETER(R20 = 3, ASYM = 3, XNEAR = HALF)
      PARAMETER(LIMIT = 20000, NDROP = 8, JMAX = 50)

#include "c309prec.inc"

      PARAMETER(HPI  = 1.57079 63267 94896 619D0)
      PARAMETER(TLOG = 0.69314 71805 59945 309D0)

      DIMENSION FC(0:*),GC(0:*),FCP(0:*),GCP(0:*),SIG(0:*)
      DIMENSION XRCF(JMAX,4)
      INTEGER NPQ(2),KAS(2)

      NINTC(W)=NINT(REAL(W))
      ABSC(W)=ABS(REAL(W))+ABS(AIMAG(W))
      NPINT(W,ACCB)=ABSC(NINTC(W)-W) .LT. ACCB .AND. REAL(W) .LT. HALF
C
      MODE=MOD(ABS(MODE1),10)
      IFCP=MOD(MODE,2) .EQ. 1
      LPR=IPR .NE. 0
      IFAIL=-2
      N11=0
      NFP=0
      KAS(1)=0
      KAS(2)=0
      NPQ(1)=0
      NPQ(2)=0
      N20=0

      ACCUR=MAX(ACCU,50*ACC8)
      ACCT=HALF*ACCUR
      ACCH=SQRT(ACCUR)
      ACCB=SQRT(ACCH)
      RERR=ACCT
      FPLMAX=LOG(FPMAX)
      FPLMIN=LOG(FPMIN)
C
      CIK=ONE
      IF(KFN .GE. 3) CIK=CI*SIGN(ONE,FPMIN-AIMAG(ZZ))
      X=ZZ*CIK
      ETA=ETA1
      IF(KFN .GT. 0) ETA=ZERO
      ETANE0=ABSC(ETA) .GT. ACC8
      ETAI=ETA*CI
      DELL=ZERO
      IF(KFN .GE. 2) DELL=HALF
      ZM1=ZLMIN-DELL
      SCALE=ZERO
      IF(MODE1 .LT. 0) SCALE=AIMAG(X)
C
      M1=1
      L1=M1+NL
      RLEL=ABS(AIMAG(ETA))+ABS(AIMAG(ZM1)) .LT. ACC8
      ABSX=ABS(X)
      AXIAL=RLEL .AND. ABS(AIMAG(X)) .LT. ACC8*ABSX
      IF(MODE .LE. 2 .AND. ABSX .LT. FPMIN) GO TO 310
      XI=ONE/X
      XLOG=LOG(X)

C       log with cut along the negative real axis, see also OMEGA

      ID=1
      DONEM=.FALSE.
      UNSTAB=.FALSE.
      LF=M1
      IFAIL=-1
   10 ZLM=ZM1+(LF-M1)
      ZLL=ZM1+(L1-M1)
C
C ***       ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels
C
      Z11=ZLL
      IF(ID .LT. 0) Z11=ZLM
      P11=CI*SIGN(ONE,ACC8-AIMAG(ETA))
      LAST=L1
C
C ***       Find phase shifts and Gamow factor at lambda = ZLL
C
      PK=ZLL+ONE
      AA=PK-ETAI
      AB=PK+ETAI
      BB=TWO*PK
      ZLNEG=NPINT(BB,ACCB)
      CLGAA=CLOGAM(AA)
      CLGAB=CLGAA
      IF(ETANE0 .AND. .NOT.RLEL) CLGAB=CLOGAM(AB)
      IF(ETANE0 .AND. RLEL) CLGAB=CONJG(CLGAA)
      SIGMA=(CLGAA-CLGAB)*CIH
      IF(KFN .EQ. 0) SIG(L1)=SIGMA
      IF(.NOT.ZLNEG) CLL=ZLL*TLOG-HPI*ETA-CLOGAM(BB)+(CLGAA+CLGAB)*HALF
      THETA=X-ETA*(XLOG+TLOG)-ZLL*HPI+SIGMA
C
      TA=(AIMAG(AA)**2+AIMAG(AB)**2+ABS(REAL(AA))+ABS(REAL(AB)))*HALF
      IF(ID .GT. 0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GO TO 20
C
C ***         use CF1 instead of CF1A, if predicted to converge faster,
C                 (otherwise using CF1A as it treats negative lambda &
C                  recurrence-unstable cases properly)
C
      RK=SIGN(ONE,REAL(X)+ACC8)
      P=THETA
      IF(RK .LT. 0) P=-X+ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA
      XRCF(1,1)=PK
      F=RK*C309R6(X*RK,ETA*RK,ZLL,P,ACCT,JMAX,NFP,FEST,ERR,FPMAX,XRCF,
     1                                      XRCF(1,3),XRCF(1,4))
      FESL=LOG(FEST)+ABS(AIMAG(X))
      NFP=-NFP
      IF(NFP .LT. 0 .OR. UNSTAB .AND. ERR .LT. ACCB) GO TO 40
      IF(.NOT.ZLNEG .OR. UNSTAB .AND. ERR .GT. ACCB) GO TO 20
      IF(LPR) WRITE(6,1060) '-L',ERR
      IF(ERR.GT.ACCB) GO TO 280
      GO TO 40
C
C ***    evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)
C
   20 IF(AXIAL) THEN
       SX1=X
       SETA=ETA
       SZLL=ZLL
       F=C309R4(SX1,SETA,SZLL,ACC8,SF,RK,ETANE0,LIMIT,ERR,NFP,
     1          FPMIN,FPMAX,LPR)
       FCL=SF
       TPK1=RK
      ELSE
       F=C309R5(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     1          FPMIN,FPMAX,LPR)
      END IF
      IF(ERR .GT. ONE) THEN
       IFAIL=-1
       GO TO 290
      END IF
C
C ***  Make a simple check for CF1 being badly unstable:
C
      IF(ID .GE. 0) THEN
       UNSTAB=REAL((ONE-ETA*XI)*CI*AIMAG(THETA)/F) .GT. ZERO
     1  .AND. .NOT.AXIAL .AND. ABS(AIMAG(THETA)) .GT. -LOG(ACC8)*HALF
     2  .AND. ABSC(ETA)+ABSC(ZLL) .LT. ABSC(X)
       IF(UNSTAB) THEN
        ID=-1
        LF=L1
        L1=M1
        RERR=ACCT
        GO TO 10
       END IF
      END IF
C
C *** compare accumulated phase FCL with asymptotic phase for G(k+1) :
C     to determine estimate of F(ZLL) (with correct sign) to start recur
C
      W=X*X*(HALF/TPK1+ONE/TPK1**2)+ETA*(ETA-TWO*X)/TPK1
      FESL=(ZLL+ONE)*XLOG+CLL-W-LOG(FCL)
   40 FESL=FESL-ABS(SCALE)
      RK=MAX(REAL(FESL),FPLMIN*HALF)
      FESL=CMPLX(MIN(RK,FPLMAX*HALF),AIMAG(FESL))
      FEST=EXP(FESL)
C
      RERR=MAX(RERR,ERR,ACC8*ABS(REAL(THETA)))
C
      FCL=FEST
      FPL=FCL*F
      IF(IFCP) FCP(L1)=FPL
      FC(L1)=FCL
C
C *** downward recurrence to lambda = ZLM. array GC,if present,stores RL
C
      I=MAX(-ID,0)
      ZL=ZLL+I
      MONO=0
      OFF=ABS(FCL)
      TA=ABSC(SIGMA)

C
C     CORRESPONDS TO   DO 70 L = L1-ID,LF,-ID
C
      L=L1-ID
      LC70=(L1-LF)/ID
   70 IF(LC70 .LE. 0) GO TO 80
      IF(ETANE0) THEN
       IF(RLEL) THEN
        DSIG=ATAN2(REAL(ETA),REAL(ZL))
        RL=SQRT(REAL(ZL)**2+REAL(ETA)**2)
       ELSE
        AA=ZL-ETAI
        BB=ZL+ETAI
        IF(ABSC(AA) .LT. ACCH .OR. ABSC(BB) .LT. ACCH) GOTO 50
        DSIG=(LOG(AA)-LOG(BB))*CIH
        RL=AA*EXP(CI*DSIG)
       END IF
       IF(ABSC(SIGMA) .LT. TA*HALF) THEN

C               re-calculate SIGMA because of accumulating roundoffs:

        SL=(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CIH
        RL=(ZL-ETAI)*EXP(CI*ID*(SIGMA-SL))
        SIGMA=SL
        TA=ZERO
       ELSE
        SIGMA=SIGMA-DSIG*ID
       END IF
       TA=MAX(TA,ABSC(SIGMA))
       SL=ETA+ZL*ZL*XI
       PL=ZERO
       IF(ABSC(ZL) .GT. ACCH) PL=(SL*SL-RL*RL)/ZL
       FCL1=(FCL*SL+ID*ZL*FPL)/RL
       SF=ABS(FCL1)
       IF(SF .GT. FPMAX) GO TO 350
       FPL=(FPL*SL+ID*PL*FCL)/RL
       IF(MODE .LE. 1) GCP(L+ID)=PL*ID
      ELSE

C                      ETA = 0, including Bessels.  NB RL==SL

       RL=ZL*XI
       FCL1=FCL*RL+FPL*ID
       SF=ABS(FCL1)
       IF(SF .GT. FPMAX) GO TO 350
       FPL=(FCL1*RL-FCL)*ID
      END IF
      MONO=MONO+1
      IF(SF .GE. OFF) MONO=0
      FCL=FCL1
      OFF=SF
      FC(L)=FCL
      IF(IFCP) FCP(L)=FPL
      IF(KFN .EQ. 0) SIG(L)=SIGMA
      IF(MODE .LE. 2) GC(L+ID)=RL
      ZL=ZL-ID
      IF(MONO .LT. NDROP .OR. AXIAL .OR.
     1             REAL(ZLM)*ID .GT. -NDROP .AND. .NOT.ETANE0) THEN
       L=L-ID
       LC70=LC70-1
       GO TO 70
      END IF
      UNSTAB=.TRUE.
C
C ***    take action if cannot or should not recur below this ZL:

   50 ZLM=ZL
      LF=L
      IF(ID .LT. 0) GO TO 380
      IF(.NOT.UNSTAB) LF=L+1
      IF(L+MONO .LT. L1-2 .OR. ID .LT. 0 .OR. .NOT.UNSTAB) GO TO 80

C             otherwise, all L values (for stability) should be done
C                        in the reverse direction:

      ID=-1
      LF=L1
      L1=M1
      RERR=ACCT
      GO TO 10

   80 IF(FCL .EQ. ZERO) FCL=ACC8
      F=FPL/FCL
C
C *** Check, if second time around, that the 'f' values agree!
C
      IF(ID .GT. 0) FIRST=F
      IF(DONEM) RERR=MAX(RERR,ABSC(F-FIRST)/ABSC(F))
      IF(DONEM) GO TO 90
C
      NOCF2=.FALSE.
      THETAM=X-ETA*(XLOG+TLOG)-ZLM*HPI+SIGMA
C
C *** on left x-plane, determine OMEGA by requiring cut on -x axis
C     on right x-plane, choose OMEGA (using estimate based on THETAM)
C       so H(omega) is smaller and recurs upwards accurately.
C     (x-plane boundary is shifted to give CF2(LH) a chance to converge)
C
      OMEGA=SIGN(ONE,AIMAG(X)+ACC8)
      IF(REAL(X) .GE. XNEAR) OMEGA=SIGN(ONE,AIMAG(THETAM)+ACC8)
      SFSH=EXP(OMEGA*SCALE-ABS(SCALE))
      OFF=EXP(MIN(TWO*MAX(ABS(AIMAG(X)),ABS(AIMAG(THETAM)),
     1                         ABS(AIMAG(ZLM))*3),FPLMAX))
      EPS=MAX(ACC8,ACCT*HALF/OFF)
C
C ***    Try first estimated omega, then its opposite,
C        to find the H(omega) linearly independent of F
C        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega)
C
   90 DO 100 L = 1,2
      LH=1
      IF(OMEGA .LT. ZERO) LH=2
      PM=CI*OMEGA
      ETAP=ETA*PM
      IF(DONEM) GO TO 130
      PQ1=ZERO
      PACCQ=ONE
      KASE=0
C
C ***            Check for small X, i.e. whether to avoid CF2 :
C
      IF(MODE .GE. 3 .AND. ABSX .LT. ONE ) GO TO 190
      IF(MODE .LT. 3 .AND. (NOCF2 .OR. ABSX .LT. XNEAR .AND.
     1   ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM) .LT. 4)) THEN
       KASE=5
       GO TO 120
      END IF
C
C ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM
C
      PQ1=C309R1(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH,
     1             LPR,ACCUR,DELL)
      ERR=ERR*MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8))
C
C *** check if impossible to get F-PQ accurately because of cancellation
C                original guess for OMEGA (based on THETAM) was wrong
C                Use KASE 5 or 6 if necessary if Re(X) < XNEAR
      IF(ERR .LT. ACCH) GO TO 110
      NOCF2=REAL(X) .LT. XNEAR .AND. ABS(AIMAG(X)) .LT. -LOG(ACC8)
  100 OMEGA=-OMEGA
      IF(UNSTAB) GO TO 360
      IF(LPR .AND. REAL(X) .LT. -XNEAR) WRITE(6,1060) '-X',ERR
  110 RERR=MAX(RERR,ERR)
C
C ***  establish case of calculation required for irregular solution
C
  120 IF(KASE .GE. 5) GO TO 130
      IF(REAL(X) .GT. XNEAR) THEN

C          estimate errors if KASE 2 or 3 were to be used:

       PACCQ=EPS*OFF*ABSC(PQ1)/MAX(ABS(AIMAG(PQ1)),ACC8)
      END IF
      IF(PACCQ .LT. ACCUR) THEN
       KASE=2
       IF(AXIAL) KASE=3
      ELSE
       KASE=1
       IF(NPQ(1)*R20 .LT. JMAX) KASE=4
C             i.e. change to kase=4 if the 2F0 predicted to converge
      END IF
  130 GO TO (190,140,150,170,190,190), ABS(KASE)
  140 IF(.NOT.DONEM)
C
C ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2)
C
     1  PQ2=C309R1(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH,
     2             LPR,ACCUR,DELL)
      P=(PQ2+PQ1)*HALF
      Q=(PQ2-PQ1)*HALF*PM
      GO TO 160
  150 P=REAL(PQ1)
      Q=AIMAG(PQ1)
C
C ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*
C
      PQ2=CONJG(PQ1)
C
C *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL
C
  160 W=(PQ1-F)*(PQ2-F)
      SF=EXP(-ABS(SCALE))
      FCM=SQRT(Q/W)*SF

C                  any SQRT given here is corrected by
C                  using sign for FCM nearest to phase of FCL

      IF(REAL(FCM/FCL) .LT. ZERO) FCM=-FCM
      GAM=(F-P)/Q
      TA=ABSC(GAM+PM)
      PACCQ=EPS*MAX(TA,ONE/TA)
      HCL=FCM*(GAM+PM)*(SFSH/(SF*SF))

C                            Consider a KASE = 1 Calculation

      IF(PACCQ .GT. ACCUR .AND. KASE .GT. 0) THEN
       F11V=C309R2(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
       IF(ERR .LT. PACCQ) GO TO 200
      END IF
      RERR=MAX(RERR,PACCQ)
      GO TO 230
C
C *** Arrive here if KASE = 4
C     to evaluate the exponentially decreasing H(LH) directly.
C
  170 IF(DONEM) GO TO 180
      AA=ETAP-ZLM
      BB=ETAP+ZLM+ONE
      F20V=C309R3(AA,BB,-HALF*PM*XI,ACCT,JMAX,ERR,FPMAX,N20,XRCF)
      IF(N20 .LE. 0) GO TO 190
      RERR=MAX(RERR,ERR)
      HCL=FPMIN
      IF(ABS(REAL(PM*THETAM)+OMEGA*SCALE) .GT. FPLMAX) GO TO 330
  180 HCL=F20V*EXP(PM*THETAM+OMEGA*SCALE)
      FCM=SFSH/((F-PQ1)*HCL)
      GO TO 230
C
C *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)
C
C           for small values of X, calculate F(X,SL) directly from 1F1
C               using REAL*16 arithmetic if possible.
C           where Z11 = ZLL if ID>0, or = ZLM if ID<0
C
  190 F11V=C309R2(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
  200 IF(N11 .LT. 0) THEN

C                               F11 failed from BB = negative integer

       IF(LPR) WRITE(6,1060) '-L',ONE
       IFAIL=-1
       GO TO 290
      END IF

C                      Consider a KASE 2 or 3 calculation :

      IF(ERR .GT. PACCQ .AND. PACCQ .LT. ACCB) THEN
       KASE=-2
       IF(AXIAL) KASE=-3
       GO TO 130
      END IF
      RERR=MAX(RERR,ERR)
      IF(ERR .GT. FPMAX) GO TO 370
      IF(ID .LT. 0) CLL=Z11*TLOG-HPI*ETA-CLOGAM(BB)
     1                        +CLOGAM(Z11+ONE+P11*ETA)-P11*SIGMA
      EK=(Z11+ONE)*XLOG-P11*X+CLL-ABS(SCALE)
      IF(ID .GT. 0) EK=EK-FESL+LOG(FCL)
      IF(REAL(EK) .GT. FPLMAX) GO TO 350
      IF(REAL(EK) .LT. FPLMIN) GO TO 340
      FCM=F11V*EXP(EK)
      IF(KASE .GE. 5) THEN
       IF(ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)) .LT. ACCH) KASE=6
C
C ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so
C ***      use an expansion for irregular soln from origin :
C
       SL=ZLM
       ZLNEG=REAL(ZLM) .LT. -ONE+ACCB
       IF(KASE .EQ. 5 .OR. ZLNEG) SL=-ZLM-ONE
       PK=SL+ONE
       AA=PK-ETAP
       AB=PK+ETAP
       BB=TWO*PK
       CLGAA=CLOGAM(AA)
       CLGAB=CLGAA
       IF(ETANE0) CLGAB=CLOGAM(AB)
       CLGBB=CLOGAM(BB)
       IF(KASE .EQ. 6 .AND. .NOT.ZLNEG) THEN
        IF(NPINT(AA,ACCUR)) CLGAA=CLGAB-TWO*PM*SIGMA
        IF(NPINT(AB,ACCUR)) CLGAB=CLGAA+TWO*PM*SIGMA
       END IF
       CLL=SL*TLOG-HPI*ETA-CLGBB+(CLGAA+CLGAB)*HALF
       DSIG=(CLGAA-CLGAB)*PM*HALF
       IF(KASE .EQ. 6) P11=-PM
       EK=PK*XLOG-P11*X+CLL-ABS(SCALE)
       SF=EXP(-ABS(SCALE))
       CHI=ZERO
       IF(.NOT.(KASE .EQ. 5 .OR. ZLNEG)) GO TO 210
C
C ***  Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)
C
C      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2
C
       CHI=SIGMA-DSIG-(ZLM-SL)*HPI
       F11V=C309R2(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),
     1             FPMAX,ACC8,ACC16)
       RERR=MAX(RERR,ERR)
       IF(KASE .EQ. 6) GO TO 210
       FESL=F11V*EXP(EK)
       FCL1=EXP(PM*CHI)*FCM
       HCL=FCL1-FESL
       RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL))
       HCL=HCL/SIN(CHI)*(SFSH/(SF*SF))
       GO TO 220
C
C *** Use the logarithmic expansion for the irregular solution (KASE 6)
C        for the case that BB is integral so sin(CHI) would be zero.
C
  210  RL=BB-ONE
       N=NINTC(RL)
       ZLOG=XLOG+TLOG-PM*HPI
       CHI=CHI+PM*THETAM+OMEGA*SCALE+AB*ZLOG
       AA=ONE-AA
       IF(NPINT(AA,ACCUR)) THEN
        HCL=ZERO
       ELSE
        IF(ID .GT. 0 .AND. .NOT.ZLNEG) F11V=FCM*EXP(-EK)
        HCL=EXP(CHI-CLGBB-CLOGAM(AA))*(-1)**(N+1)*(F11V*ZLOG+
     1   C309R2(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16))
        RERR=MAX(RERR,ERR)
       END IF
       IF(N .GT. 0) THEN
        EK=CHI+CLOGAM(RL)-CLGAB-RL*ZLOG
        DF=C309R2(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16)
        HCL=HCL+EXP(EK)*DF
       END IF
       RERR=MAX(RERR,TWO*ABS(BB-NINTC(BB)))
  220  PQ1=F-SFSH/(FCM*HCL)
      ELSE
       IF(MODE .LE. 2) HCL=SFSH/((F-PQ1)*FCM)
       KASE=1
      END IF
C
C ***  Now have absolute normalisations for Coulomb Functions
C          FCM & HCL  at lambda = ZLM
C      so determine linear transformations for Functions required :
C
  230 IH=ABS(MODE1)/10
      IF(KFN .EQ. 3) IH=(3-AIMAG(CIK))/2+HALF
      P11=ONE
      IF(IH .EQ. 1) P11=CI
      IF(IH .EQ. 2) P11=-CI
      DF=-PM
      IF(IH .GE. 1) DF=-PM+P11
      IF(ABSC(DF) .LT. ACCH) DF=ZERO
C
C *** Normalisations for spherical or cylindrical Bessel functions
C
      IF(KFN .LE. 0) THEN
       ALFA=ZERO
       BETA=ONE
      ELSE IF(KFN .EQ. 1) THEN
       ALFA=XI
       BETA=XI
      ELSE
       ALFA=XI*HALF
       BETA=SQRT(XI/HPI)
       IF(REAL(BETA) .LT. ZERO) BETA=-BETA
      END IF
      AA=ONE
      IF(KFN .GT. 0) AA=-P11*BETA

C                Calculate rescaling factors for I & K output

      IF(KFN .GE. 3) THEN
       P=EXP((ZLM+DELL)*HPI*CIK)
       AA=BETA*HPI*P
       BETA=BETA/P
       Q=CIK*ID
      END IF

C                  Calculate rescaling factors for GC output

      IF(IH .EQ. 0) THEN
       TA=ABS(SCALE)+AIMAG(PM)*SCALE
       RK=ZERO
       IF(TA .LT. FPLMAX) RK=EXP(-TA)
      ELSE
       TA=ABS(SCALE)+AIMAG(P11)*SCALE
       IF(ABSC(DF) .GT. ACCH .AND. TA .GT. FPLMAX) GO TO 320
       IF(ABSC(DF) .GT. ACCH) DF=DF*EXP(TA)
       SF=TWO*(LH-IH)*SCALE
       RK=ZERO
       IF(SF .GT. FPLMAX) GO TO 320
       IF(SF .GT. FPLMIN) RK=EXP(SF)
      END IF
      KAS((3-ID)/2)=KASE
      W=FCM/FCL
      IF(LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GO TO 340
      IF(MODE .GE. 3) GO TO 240
      IF(LPR .AND. ABSC(F-PQ1) .LT. ACCH*ABSC(F))
     1                             WRITE(6,1020) LH,ZLM+DELL
      HPL=HCL*PQ1
      IF(ABSC(HPL) .LT. FPMIN .OR. ABSC(HCL) .LT. FPMIN) GO TO 330
C
C *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)
C *** renormalise FC,FCP at each lambda
C ***    ZL   = ZLM - MIN(ID,0) here
C
  240 DO 270 L = LF,L1,ID
      FCL=W*FC(L)
      IF(ABSC(FCL) .LT. FPMIN) GO TO 340
      IF(IFCP) FPL=W*FCP(L)
      FC(L)=BETA*FCL
      IF(IFCP) FCP(L)=BETA*(FPL-ALFA*FCL)*CIK
      FC(L)=C309R8(FC(L),ACCUR)
      IF(IFCP) FCP(L)=C309R8(FCP(L),ACCUR)
      IF(MODE .GE. 3) GO TO 260
      IF(L .EQ. LF) GO TO 250
      ZL=ZL+ID
      ZID=ZL*ID
      RL=GC(L)
      IF(ETANE0) THEN
       SL=ETA+ZL*ZL*XI
       IF(MODE .EQ. 1) THEN
        PL=GCP(L)
       ELSE
        PL=ZERO
        IF(ABSC(ZL) .GT. ACCH) PL=(SL*SL-RL*RL)/ZID
       END IF
       HCL1=(SL*HCL-ZID*HPL)/RL
       HPL=(SL*HPL-PL*HCL)/RL
      ELSE
       HCL1=RL*HCL-HPL*ID
       HPL=(HCL-RL*HCL1)*ID
      END IF
      HCL=HCL1
      IF(ABSC(HCL) .GT. FPMAX) GO TO 320
  250 GC(L)=AA*(RK*HCL+DF*FCL)
      IF(MODE .EQ. 1) GCP(L)=(AA*(RK*HPL+DF*FPL)-ALFA*GC(L))*CIK
      GC(L)=C309R8(GC(L),ACCUR)
      IF(MODE .EQ. 1) GCP(L)=C309R8(GCP(L),ACCUR)
      IF(KFN .GE. 3) AA=AA*Q
  260 IF(KFN .GE. 3) BETA=-BETA*Q
  270 LAST=MIN(LAST,(L1-L)*ID)
      GO TO 280
C
C *** Come here after all soft errors to determine how many L values ok
C
  310 IF(LPR) WRITE(6,1000) ZZ
      GO TO 999
  320 IF(LPR) WRITE(6,1010) ZL+DELL,'IR',HCL,'>',FPMAX
      GO TO 280
  330 IF(LPR) WRITE(6,1010) ZL+DELL,'IR',HCL,'<',FPMIN
      GO TO 280
  340 IF(LPR) WRITE(6,1010) ZL+DELL,'  ',FCL,'<',FPMIN
      GO TO 280
  350 IF(LPR) WRITE(6,1010) ZL+DELL,'  ',FCL,'>',FPMAX
      GO TO 280
  360 IF(LPR) WRITE(6,1030) ZL+DELL
      GO TO 280
  370 IF(LPR) WRITE(6,1040) Z11,I
      IFAIL=-1
      GO TO 290
  380 IF(LPR) WRITE(6,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL
      IFAIL=-1
      GO TO 290
  280 IF(ID .GT. 0 .OR. LAST .EQ. 0) IFAIL=LAST
      IF(ID .LT. 0 .AND. LAST .NE. 0) IFAIL=-3
C
C *** Come here after ALL errors for this L range (ZLM,ZLL)
C
C *** so on first block, 'F' started decreasing monotonically,
C                        or hit bound states for low ZL.
C     thus redo M1 to LF-1 in reverse direction, i.e. do
C      CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)
C
  290 IF(ID .GT. 0 .AND. LF .NE. M1) THEN
       ID=-1
       IF(.NOT.UNSTAB) LF=LF-1
       DONEM=UNSTAB
       LF=MIN(LF,L1)
       L1=M1
       GO TO 10
      END IF
      IF(IFAIL .LT. 0) GO TO 999
      IF(LPR .AND. RERR .GT. ACCB) WRITE(6,1070) RERR
      IF(RERR .GT. 0.1D0) IFAIL=-4
  999 DO 998 L = 1,NL+1
      FC(L-1)=FC(L)
      GC(L-1)=GC(L)
      FCP(L-1)=FCP(L)
      GCP(L-1)=GCP(L)
  998 SIG(L-1)=SIG(L)
      RETURN
C
 1000 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',
     2 1P2D10.2,' ABS(X) TOO SMALL')
 1010 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P2E10.1,') ',
     2 A1,E10.1)
 1020 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''H(',I1,
     2 ')'' IS LOST AT ZL = ',2F7.2/1X,'*****',22X,'(EG. COULOMB ',
     3 'EIGENSTATE OR CF1 UNSTABLE)')
 1030 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     2 '(ETA & L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE AT L = ',2F8.2)
 1040 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'OVERFLOW IN 1F1 SERIES AT ZL = ',2F8.3,' AT TERM',I5)
 1050 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR OR MULTIPLE',
     2 ' INSTABILITIES PRESENT'/
     3   1X,'*****',22X,'TRY CALLING TWICE, 1ST FOR ZL FROM',2F8.3,
     4 ' TO',2F8.3,' (INCL)'/1X,'*****',41X,'2ND FOR ZL FROM',2F8.3,
     5 ' TO',2F8.3)
 1060 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'WARNING: AS ''',A2,''' REFLECTION RULES NOT USED ',
     2 'ERRORS CAN BE UP TO',1PD12.2)
 1070 FORMAT(1X,'***** CERN C309 CCLBES ... ',
     1 'WARNING: OVERALL ROUNDOFF ERROR APPROXIMATELY',1PE11.1)
      END
#endif

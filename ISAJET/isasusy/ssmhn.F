#include "isajet/pilot.h"
      SUBROUTINE SSMHN(MHLNEG)
C-----------------------------------------------------------------------
C
C          Calculate HL, HH masses and ALFAH 
C          (scalar Higgs mixing angle) using radiative
C          corrections calculated by M. Bisset
C          and save results in /SSPAR/.
C
C          Both top and bottom couplings are now 
C          included.  Non-degenerate mixed squark
C          masses and A-terms are also included.
C          The D-terms from the squark mass matrix
C          (terms prop. to g**2 * Yukawa coupling)
C          are included as an option: 
C                 INRAD = 1 ==> D-TERMS ON
C                 INRAD = 2 ==> D-TERMS OFF    .
C
C         10/18/93 D-terms are now turned on.
C                     INRAD = 1 
C
C         There is an arbitrary mass scale that must
C         chosen to avoid dimensionful logarithms.
C         The choice does not matter if D-terms are
C         not included, but it does matter if D-terms
C         are included. 
C     
C         Arbitrary mass scale updated to 
C              QQQ = HIGFRZ = SQRT(AMTLSS*AMTRSS)
C         with running masses to include dominant 2-loop 
C         effects. 12/10/96 H. Baer
C
C         It is assumed that the A-terms are real.
C         Complex A-terms are allowed 
C         (unless RTT=0 or RBB=0 --see below) in 
C         this subroutine, but the imaginary parts
C         are now set to zero. 
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
C
      REAL PI,PI2,SR2,G2,GP2,GGP,GG1,GG2
      REAL TANB,COTB,COSB,SINB,BE
      REAL SINB2,COSB2,COS2B,SIN2B
      REAL V2,VP2,V,VP,VVP,VPVM,VVPP
      REAL MT2,MB2,FT2,FB2,MW2,ZAP,QQQ2
      REAL EP,EP2,RR,MHP2
      REAL ATI,ABI,ATR,ABR,AT2,AB2
      REAL TLRM,BLRM,TLRP,BLRP
      REAL MST1SQ,MST2SQ,MSB1SQ,MSB2SQ

      REAL RTT,TTT1,TEMPT,TM1BT
      REAL TEMPS,T1RD,T2RD,T1RPD,T2RPD
      REAL CT1,A1,A2,T1RR,T2RR
      REAL CT5,A5,A6,T1RPRP,T2RPRP
      REAL A9,T1RRP,T2RRP
      REAL TEMPSQ,DT1,DT2,VRRT,VRPRPT,VRRPT
      REAL ALPHAT,LAT
C
      REAL RBB,BBB1,TEMPB,TM1BB
      REAL B1RD,B2RD,B1RPD,B2RPD
      REAL CB3,A3,A4,B1RR,B2RR
      REAL CB7,A7,A8,B1RPRP,B2RPRP
      REAL A10,B1RRP,B2RRP
      REAL DB1,DB2,VRRB,VRPRPB,VRRPB
      REAL ALPHAB,LAB
C
      REAL DVRR,DVRPRP,DVRRP,TEMPH
      REAL MHL2,MHH2,TRACEM,TPAL,TANAH
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS,HIGFRZ
      DOUBLE PRECISION SSMQCD
      INTEGER INRAD,MHLNEG
C
      MHLNEG=0
      PI=4.*ATAN(1.)
      PI2 = PI**2
      SR2=SQRT(2.)
      G2=4.*PI*ALFAEM/SN2THW
      GP2=G2*SN2THW/(1.-SN2THW)
      HIGFRZ=SQRT(AMTLSS*AMTRSS)
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(HIGFRZ))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(HIGFRZ))
      MT2=MTQ**2
      MB2=MBQ**2
      MW2=AMW**2
      EP=TWOM1
      EP2=EP**2
      RR=RV2V1
      MHP2=AMHA**2
      TANB=1.0/RR
      COTB=RR
      BE=ATAN(1./RV2V1)
      SINB=SIN(BE)
      COSB=COS(BE)
      SINB2=SINB**2
      COSB2=COSB**2
      SIN2B=SIN(2.*BE)
      COS2B=COS(2.*BE)
      V2=2.0*MW2*SINB2/G2
      VP2=2.0*MW2*COSB2/G2
      V=SQRT(V2)
      VP=SQRT(VP2)
      VVP=SQRT(V2*VP2)
      VPVM=VP2-V2
      GGP=G2+GP2
      GG1=G2-5.0*GP2/3.0
      GG2=G2-GP2/3.0
      VVPP=2.0*AMZ**2/GGP
      FT2=MT2/V2
      FB2=MB2/VP2
C
      TLRM=AMTLSS**2-AMTRSS**2
      BLRM=AMBLSS**2-AMBRSS**2 
      TLRP=AMTLSS**2+AMTRSS**2
      BLRP=AMBLSS**2+AMBRSS**2 
C
C          Higgs mass matrix
C
C     (AAT and AAB are also assumed to be real)
C          
      ATR=AAT
      ABR=AAB
      ATI=0.0
      ABI=0.0
      AT2=ATR**2+ATI**2
      AB2=ABR**2+ABI**2
C
      MST1SQ=AMT1SS**2
      MST2SQ=AMT2SS**2
      MSB1SQ=AMB1SS**2
      MSB2SQ=AMB2SS**2
      INRAD=1
      QQQ2=HIGFRZ**2
C
      ZAP = 1.0
C
C                  STOP TERMS
C          
      RTT=(TLRM+VPVM*ZAP*GG1/4.0)**2
     $      +4.0*MT2*(EP*COTB+ATR)**2+4.0*MT2*ATI**2
      RTT=SQRT(RTT)
C
C     calculate 2M1*B term
C
      TTT1=0.5*TLRP+MT2+VPVM*ZAP*GGP/8.0
      IF(RTT.NE.0.0) THEN
        TEMPT=4.0*EP*FT2*VVP*ATI**2/(RTT**2)
        TM1BT=-2.0*FT2*(TEMPT+ATR)*TTT1
     $               *LOG(MST2SQ/MST1SQ)/RTT
        TM1BT=TM1BT-FT2*ATR
     $               *LOG(MST1SQ*MST2SQ/QQQ2/QQQ2)
        TM1BT=TM1BT+FT2*(2.0*TEMPT-ATR)
        TM1BT=3.0*EP*TM1BT/32.0/PI2
C
C        calculate first derivatives w.r.t H_R
C           divided by sqrt(2) * v
C        
         TEMPS=-ZAP*GG1*(TLRM+ZAP*GG1*VPVM/4.0)/2.0 
         TEMPS=TEMPS+4.0*FT2*(AT2+EP*COTB*ATR)
         TEMPS=TEMPS/RTT/4.0 
         T1RD=FT2-ZAP*GGP/8.0-TEMPS
         T2RD=FT2-ZAP*GGP/8.0+TEMPS
C
C        calculate first derivatives w.r.t H_R'
C           divided by sqrt(2) * v'
C        
         TEMPS=ZAP*GG1*(TLRM+ZAP*GG1*VPVM/4.0)/2.0
         TEMPS=TEMPS+4.0*FT2*EP*(EP+TANB*ATR)
         TEMPS=TEMPS/RTT/4.0
         T1RPD=ZAP*GGP/8.0-TEMPS
         T2RPD=ZAP*GGP/8.0+TEMPS
C
C        calculate second derivatives w.r.t. H_R
C
         CT1=-V*ZAP*GG1*(TLRM+ZAP*GG1*VPVM/4.0)/SR2
         CT1=CT1+4.0*SR2*FT2*V*(EP*COTB*ATR+AT2)
         A1=-CT1**2/(RTT**3)/8.0
         A2=-ZAP*GG1*(TLRM+ZAP*GG1*VPVM/4.0)/2.0
         A2=A2+V2*ZAP*GG1**2/4.0+4.0*FT2*AT2
         A2=A2/RTT/4.0
         T1RR=FT2-ZAP*GGP/8.0-A1-A2
         T2RR=FT2-ZAP*GGP/8.0+A1+A2
C
C        calculate second derivatives w.r.t. H_R'
C
         CT5=VP*ZAP*GG1*(TLRM+ZAP*GG1*VPVM/4.0)/SR2
         CT5=CT5+4.0*SR2*FT2*VP*EP*(EP+TANB*ATR)
         A5=-CT5**2/(RTT**3)/8.0
         A6=ZAP*GG1*(TLRM+ZAP*GG1*VPVM/4.0)/2.0 
         A6=A6+VP2*ZAP*GG1**2/4.0+4.0*FT2*EP2
         A6=A6/RTT/4.0
         T1RPRP=ZAP*GGP/8.0-A5-A6
         T2RPRP=ZAP*GGP/8.0+A5+A6
C
C        calculate second derivatives w.r.t. H_R and H_R'
C
         A9=-VVP*ZAP*(GG1**2)/4.0+4.0*FT2*EP*ATR
         A9=A9/RTT/4.0
         T1RRP=CT1*CT5/(RTT**3)/8.0-A9
         T2RRP=-CT1*CT5/(RTT**3)/8.0+A9
C
C        calculate D^2 V / D^2 H_R
C
         TEMPSQ=MST1SQ*(T1RR-T1RD)
         DT1=2.0*(2.0*V2*T1RD**2+TEMPSQ)*LOG(MST1SQ/QQQ2)
         DT1=DT1+6.0*V2*T1RD**2+TEMPSQ
         TEMPSQ=MST2SQ*(T2RR-T2RD)
         DT2=2.0*(2.0*V2*T2RD**2+TEMPSQ)*LOG(MST2SQ/QQQ2)
         DT2=DT2+6.0*V2*T2RD**2+TEMPSQ
         VRRT=DT1+DT2-8.0*FT2*MT2*LOG(MT2/QQQ2)-12.0*FT2*MT2
         VRRT=-TM1BT*COTB+3.0*VRRT/32.0/PI2
C
C        calculate D^2 V / D^2 H'_R
C
         TEMPSQ=MST1SQ*(T1RPRP-T1RPD)
         DT1=2.0*(2.0*VP2*T1RPD**2+TEMPSQ)*LOG(MST1SQ/QQQ2)
         DT1=DT1+6.0*VP2*T1RPD**2+TEMPSQ
         TEMPSQ=MST2SQ*(T2RPRP-T2RPD)
         DT2=2.0*(2.0*VP2*T2RPD**2+TEMPSQ)*LOG(MST2SQ/QQQ2)
         DT2=DT2+6.0*VP2*T2RPD**2+TEMPSQ
         VRPRPT=-TM1BT*TANB+3.0*(DT1+DT2)/32.0/PI2
C
C        calculate D^2 V / D^H_R  D^H_R'
C
         DT1=2.0*VVP*T1RD*T1RPD+MST1SQ*T1RRP
         DT1=2.0*DT1*LOG(MST1SQ/QQQ2)
         DT1=DT1+6.0*VVP*T1RD*T1RPD+MST1SQ*T1RRP
         DT2=2.0*VVP*T2RD*T2RPD+MST2SQ*T2RRP
         DT2=2.0*DT2*LOG(MST2SQ/QQQ2)
         DT2=DT2+6.0*VVP*T2RD*T2RPD+MST2SQ*T2RRP
         VRRPT=TM1BT+3.0*(DT1+DT2)/32.0/PI2
C
      ELSE IF(RTT.EQ.0.0) THEN
C
         ALPHAT=TLRP/2.0+MT2+ZAP*GGP*VPVM/8.0
         LAT=2.0*LOG(ALPHAT/QQQ2)+3.0
C
C        calculate D^2 V / D^2 H_R
C
         VRRT=V2*(GGP**2+GG1**2)/16.0-MT2*GGP
         VRRT=ZAP*VRRT*LAT+8.0*FT2*MT2*LOG(ALPHAT/MT2)
         VRRT=3.0*VRRT/32.0/PI2
C
C        calculate D^2 V / D^2 H_R'
C
         VRPRPT=ZAP*VP2*(GGP**2+GG1**2)/16.0
         VRPRPT=3.0*(VRPRPT*LAT)/32.0/PI2
C
C        calculate D^2 V / D^H_R D^H_R'
C
         VRRPT=FT2*GGP-(GGP**2+GG1**2)/8.0
         VRRPT=ZAP*VVP*VRRPT*LAT/2.0
         VRRPT=3.0*VRRPT/32.0/PI2
C
C
      ENDIF
C
C     SBOTTOM TERMS
C
      RBB=(BLRM-VPVM*ZAP*GG2/4.0)**2
     $      +4.0*MB2*(EP*TANB+ABR)**2+4.0*MB2*ABI**2
      RBB=SQRT(RBB)
C      IF(RBB.EQ.0.0.AND.ABI.NE.0.0) THEN
C        WRITE(6,*) 'RBB=0, ABI NOT 0'
C        WRITE(6,*) 'ERROR: THIS CASE NOT COVERED YET'
C        GO TO 1000
C      ENDIF
C
      IF(RBB.NE.0.0) THEN
C
C     calculate 2M1*B term
C
        BBB1=0.5*BLRP+MB2-VPVM*ZAP*GGP/8.0
        TEMPB=4.0*EP*FB2*VVP*ABI**2/(RBB**2)       
        TM1BB=-2.0*FB2*(TEMPB+ABR)*BBB1
     $               *LOG(MSB2SQ/MSB1SQ)/RBB
        TM1BB=TM1BB-FB2*ABR
     $               *LOG(MSB1SQ*MSB2SQ/QQQ2/QQQ2)
        TM1BB=TM1BB+FB2*(2.0*TEMPB-ABR)
        TM1BB=3.0*EP*TM1BB/32.0/PI2
C
C        calculate first derivatives w.r.t H_R
C           divided by sqrt(2) * v
C        
        TEMPS=ZAP*GG2*(BLRM-ZAP*GG2*VPVM/4.0)/2.0
        TEMPS=TEMPS+4.0*FB2*EP*(EP+COTB*ABR)
        TEMPS=TEMPS/RBB/4.0
        B1RD=ZAP*GGP/8.0-TEMPS
        B2RD=ZAP*GGP/8.0+TEMPS

C        calculate first derivatives w.r.t H_R'
C           divided by sqrt(2) * v'
C
        TEMPS=-ZAP*GG2*(BLRM-ZAP*GG2*VPVM/4.0)/2.0
        TEMPS=TEMPS+4.0*FB2*(AB2+EP*TANB*ABR)
        TEMPS=TEMPS/RBB/4.0
        B1RPD=FB2-ZAP*GGP/8.0-TEMPS
        B2RPD=FB2-ZAP*GGP/8.0+TEMPS
C
C        calculate second derivatives w.r.t. H_R
C
        CB3=V*ZAP*GG2*(BLRM-ZAP*GG2*VPVM/4.0)/SR2
        CB3=CB3+4.0*SR2*FB2*V*EP*(EP+COTB*ABR)
        A3=-CB3**2/(RBB**3)/8.0
        A4=ZAP*GG2*(BLRM-ZAP*GG2*VPVM/4.0)/2.0
        A4=A4+V2*ZAP*GG2**2/4.0+4.0*FB2*EP2
        A4=A4/RBB/4.0
        B1RR=ZAP*GGP/8.0-A3-A4
        B2RR=ZAP*GGP/8.0+A3+A4
C
C       calculate second derivatives w.r.t. H_R'
C
        CB7=-VP*ZAP*GG2*(BLRM-ZAP*GG2*VPVM/4.0)/SR2
        CB7=CB7+4.0*SR2*FB2*VP*(AB2+EP*TANB*ABR)
        A7=-CB7**2/(RBB**3)/8.0
        A8=-ZAP*GG2*(BLRM-ZAP*GG2*VPVM/4.0)/2.0
        A8=A8+VP2*ZAP*GG2**2/4.0+4.0*FB2*AB2
        A8=A8/RBB/4.0
        B1RPRP=FB2-ZAP*GGP/8.0-A7-A8
        B2RPRP=FB2-ZAP*GGP/8.0+A7+A8
C
C       calculate second derivatives w.r.t. H_R and H_R'
C
        A10=-VVP*ZAP*(GG2**2)/4.0+4.0*FB2*EP*ABR
        A10=A10/RBB/4.0
        B1RRP=CB3*CB7/(RBB**3)/8.0-A10
        B2RRP=-CB3*CB7/(RBB**3)/8.0+A10
C
C       calculate  D^2 V / D^2 H_R
C
        TEMPSQ=MSB1SQ*(B1RR-B1RD)
        DB1=2.0*(2.0*V2*B1RD**2+TEMPSQ)*LOG(MSB1SQ/QQQ2)
        DB1=DB1+6.0*V2*B1RD**2+TEMPSQ
        TEMPSQ=MSB2SQ*(B2RR-B2RD)
        DB2=2.0*(2.0*V2*B2RD**2+TEMPSQ)*LOG(MSB2SQ/QQQ2)
        DB2=DB2+6.0*V2*B2RD**2+TEMPSQ
        VRRB=-TM1BB*COTB+3.0*(DB1+DB2)/32.0/PI2
C
C       calculate  D^2 V / D^2 H'_R
C
        TEMPSQ=MSB1SQ*(B1RPRP-B1RPD)
        DB1=2.0*(2.0*VP2*B1RPD**2+TEMPSQ)*LOG(MSB1SQ/QQQ2)
        DB1=DB1+6.0*VP2*B1RPD**2+TEMPSQ
        TEMPSQ=MSB2SQ*(B2RPRP-B2RPD)
        DB2=2.0*(2.0*VP2*B2RPD**2+TEMPSQ)*LOG(MSB2SQ/QQQ2)
        DB2=DB2+6.0*VP2*B2RPD**2+TEMPSQ
        VRPRPB=DB1+DB2
        VRPRPB=DB1+DB2-8.0*FB2*MB2*LOG(MB2/QQQ2)-12.0*FB2*MB2
        VRPRPB=-TM1BB*TANB+3.0*VRPRPB/32.0/PI2
C
C       calculate  D^2 V / D H_R  D H'_R
C
        DB1=2.0*VVP*B1RD*B1RPD+MSB1SQ*B1RRP
        DB1=2.0*DB1*LOG(MSB1SQ/QQQ2)
        DB1=DB1+6.0*VVP*B1RD*B1RPD+MSB1SQ*B1RRP
        DB2=2.0*VVP*B2RD*B2RPD+MSB2SQ*B2RRP
        DB2=2.0*DB2*LOG(MSB2SQ/QQQ2)
        DB2=DB2+6.0*VVP*B2RD*B2RPD+MSB2SQ*B2RRP
        VRRPB=TM1BB+3.0*(DB1+DB2)/32.0/PI2
 
      ELSE IF(RBB.EQ.0.0) THEN
C
        ALPHAB=BLRP/2.0+MB2-ZAP*GGP*VPVM/8.0
        LAB=2.0*LOG(ALPHAB/QQQ2)+3.0
C
C       calculate  D^2 V / D^2 H_R
C
        VRRB=ZAP*V2*(GGP**2 + GG2**2)/16.0
        VRRB=3.0*(VRRB*LAB)/32.0/PI2
C
C       calculate  D^2 V / D^2 H_R'
C
        VRPRPB=VP2*(GGP**2+GG2**2)/16.0-MB2*GGP         
        VRPRPB=ZAP*VRPRPB*LAB+8.0*FB2*MB2*LOG(ALPHAB/MB2)
        VRPRPB=3.0*VRPRPB/32.0/PI2
C
C       calculate  D^2 V / D^H_R D^H_R'
C
        VRRPB=FB2*GGP-(GGP**2+GG2**2)/8.0
        VRRPB=ZAP*VVP*VRRPB*LAB/2.0
        VRRPB=3.0*VRRPB/32.0/PI2
C
      ENDIF
C
      DVRR=VRRT+VRRB+VP2*MHP2/VVPP + V2*GGP/2.0
      DVRPRP=VRPRPT+VRPRPB+V2*MHP2/VVPP + VP2*GGP/2.0
      DVRRP=VRRPT+VRRPB-VVP*MHP2/VVPP - VVP*GGP/2.0
C          TEMPH is always non-negative:
      TEMPH=(DVRR-DVRPRP)**2+4*DVRRP**2
      TEMPH=0.5*SQRT(TEMPH)
      MHL2=0.5*(DVRR+DVRPRP)-TEMPH
      MHH2=0.5*(DVRR+DVRPRP)+TEMPH
      IF(MHL2.LT.0.0) THEN
        MHLNEG=1
C        WRITE(LOUT,*) 'SSMHN: ERROR:  MHL**2 < 0.0 FOR PARAMETERS:'
C        WRITE(LOUT,*) 'MHP =', AMHA, 'TANB =', 1.0/RR
C        WRITE(LOUT,*) 'MSTL=', AMTLSS, 'MSBL=', AMBLSS 
C        WRITE(LOUT,*) 'MSTR=', AMTRSS, 'MSBR=', AMBRSS
C        WRITE(LOUT,*) 'AT=', AAT, 'AB=', AAB
C        WRITE(LOUT,*) 'MU=-2M1=', -EP
C        WRITE(LOUT,*) 'MT=', AMTP, 'MB=', AMBT
C        WRITE(LOUT,*) 'D-TERMS? 1=YES 2=NO :', INRAD
C        WRITE(LOUT,*) 'MASS SCALE (QQQ)=', SQRT(QQQ2)
        AMHH=SQRT(MHH2)
        AMHL=SQRT(ABS(MHL2))
        GO TO 1000
      ENDIF
      AMHL=SQRT(MHL2)
      AMHH=SQRT(MHH2)

C
C     Now calculate mixing angle ALFAH
C
      TRACEM=DVRR-DVRPRP
      TPAL=TRACEM**2 + 4.0*DVRRP**2
      TANAH=TRACEM+SQRT(TPAL)
      IF(DVRRP.EQ.0.0) THEN
        WRITE(LOUT,*) 'SSMHN: OFF-DIAGONAL TERM OF SCALAR HIGGS',
     $  ' MASS MATRIX IS ZERO '
        IF(TANAH.NE.0.0) THEN
          WRITE(LOUT,*) 'SSMHN: WARNING: TAN(ALFAH) FORMULA',
     $    ' YIELDS INFINITY'
        ELSE IF(TANAH.EQ.0.0) THEN
          WRITE(LOUT,*) 'SSMHN: WARNING: TAN(ALFAH) FORMULA',
     $    ' YIELDS 0/0 '
        ENDIF
        IF(DVRR.GT.DVRPRP) THEN
          WRITE(LOUT,*) 'SSMHN: DVRR > DVRPRP ==> SET ALFAH=PI/2'
          ALFAH = PI/2.0
        ELSE IF (DVRR .LT. DVRPRP) THEN
          WRITE(LOUT,*) 'SSMHN: DVRR < DVRPRP ==> SET ALFAH=0'
          ALFAH = 0.0
        ELSE IF (DVRR .EQ. DVRPRP) THEN
          WRITE(LOUT,*) 'SSMHN: DVRR = DVRPRP ==> ALFAH INDETERMINANT'
          WRITE(LOUT,*) 'SETTING SCALAR MIXING ANGLE ALPHA=PI/4'
          ALFAH=PI/4.0
        ENDIF
        GO TO 1000
      ENDIF
      TANAH = -0.5*TANAH/DVRRP
      ALFAH = ATAN(TANAH)
C
1000  RETURN
      END

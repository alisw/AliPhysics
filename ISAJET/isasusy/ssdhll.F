#include "isajet/pilot.h"
      SUBROUTINE SSDHLL(DELHLL)
C-----------------------------------------------------------------------
C
C          Calculates radiative correction to the
C          H_h-H_l-H_l vertex.
C          calculated by M. Bisset
C
C          This subroutine calculates the 
C          radiative correction to the 
C          H_h-H_l-H_l vertex which can be 
C          important in determining the 
C           H_h --> H_l H_l partial decay width.
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
C         10/18/93 arbitrary mass scale set to H_h mass
C                       QQQ = AMHH
C
C         It is assumed that the A-terms are real.
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
      REAL V2,VP2,V,VP,VVP,VPVM,VVPP,MT,MB
      REAL MT2,MB2,FT2,FB2,FT,FB,FT4,FB4
      REAL MW2,ZAP,QQQ2,EP,EP2,RR,MHP2
      REAL ATI,ABI,ATR,ABR,AT2,AB2
      REAL MSTL2,MSTR2,MSBL2,MSBR2
      REAL TLRM,BLRM
      REAL MST1SQ,MST2SQ,MSB1SQ,MSB2SQ
      REAL RTT,RBB
C
      REAL A0,A1,A2,A1P,A2P,A3,A4
      REAL B0,B1,B2,B1P,B2P,B3,B4
      REAL MT1R,MT2R,MB1R,MB2R
      REAL MT1P,MT2P,MB1P,MB2P
      REAL MT1RR,MT2RR,MB1RR,MB2RR
      REAL MT1PP,MT2PP,MB1PP,MB2PP
      REAL MT1RP,MT2RP,MB1RP,MB2RP
      REAL MT1RRR,MT2RRR,MB1RRR,MB2RRR
      REAL MT1PRR,MT2PRR,MB1PRR,MB2PRR
      REAL MT1RPP,MT2RPP,MB1RPP,MB2RPP
      REAL MT1PPP,MT2PPP,MB1PPP,MB2PPP
C
      REAL SQVT1,SQVT2,SQVB1,SQVB2
      REAL SQVRRR,SQVPPP,SQVPRR,SQVRPP
      REAL FVRRR,FVPPP
      REAL VRRR,VPPP,VPRR,VRPP
C
      REAL ALPHAT,GGP1SQ,ALPHAB,GGP2SQ,TEMPSQ,BSQ
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS,HIGFRZ
      DOUBLE PRECISION SSMQCD
C
      REAL CA2,SA2,DVHLL
      DOUBLE PRECISION DELHLL
C
      INTEGER INRAD,ISPECT,ISPECB
C
      MW2=AMW**2
      HIGFRZ=SQRT(AMTLSS*AMTRSS)
      QQQ2=HIGFRZ**2
      INRAD=1
      ZAP=1.0
C
      PI=4.*ATAN(1.)
      PI2=PI**2
      SR2=SQRT(2.)
      G2=4.*PI*ALFAEM/SN2THW
      GP2=G2*SN2THW/(1.-SN2THW)
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(HIGFRZ))
      HIGFRZ=SQRT(AMTLSS*AMTRSS)
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(HIGFRZ))
      MT=MTQ
      MB=MBQ
      MT2=MT**2
      MB2=MB**2
      EP=TWOM1
      EP2=EP**2
      MHP2=AMHA**2
      RR=RV2V1
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
      FT=SQRT(FT2)
      FB=SQRT(FB2)
      FT4 = FT2**2
      FB4 = FB2**2
C
C      (AAT and AAB are also assumed to be real)
C
      ATR=AAT
      ABR=AAB
      ATI=0.0
      ABI=0.0
      AT2=ATR**2+ATI**2
      AB2=ABR**2+ABI**2
C
      MSTL2=AMTLSS**2
      MSTR2=AMTRSS**2
      MSBL2=AMBLSS**2
      MSBR2=AMBRSS**2
      TLRM=MSTL2-MSTR2
      BLRM=MSBL2-MSBR2 
C
C      UNFORTUNATELY, I HAVE USED MY OLD CONVENTION
C      FOR THE STOP AND SBOTTOM EIGENVALUES HERE 
C      (T1 <==> T2 OF NOTATION IN X. TATA'S AND OTHER
C      PEOPLE'S NOTATION).  SO THE NEXT FOUR LINES ARE
C      A FIX-UP UNTIL I GO THROUGH AND CHANGE THE
C      NOTATION THROUGHOUT THIS SUBROUTINE.
C
      MST2SQ=AMT1SS**2
      MST1SQ=AMT2SS**2
      MSB2SQ=AMB1SS**2
      MSB1SQ=AMB2SS**2
C
C
C      Calculation of radiative correction to 
C      the H_H-H_l-H_l vertex
C
C
C                  STOP TERMS
C
      ISPECT=0          
      RTT=(TLRM+VPVM*ZAP*GG1/4.0)**2
     $      +4.0*MT2*(EP*COTB+ATR)**2+4.0*MT2*ATI**2
C
      IF(RTT.GT.0.0) THEN
        A0=SQRT(RTT)
        A1=-V*ZAP*GG1*(TLRM+ZAP*VPVM*GG1/4.0)/SR2
        A1=A1+4.0*SR2*FT*MT*(AT2+EP*ATR*COTB)        
        A2=-ZAP*GG1*(TLRM+ZAP*VPVM*GG1/4.0)/2.0
        A2=A2 +V2*ZAP*GG1**2/4.0 +4.0*FT2*AT2
        A1P=VP*ZAP*GG1*(TLRM+ZAP*VPVM*GG1/4.0)/SR2
        A1P=A1P+4.0*SR2*FT*MT*EP*(ATR+EP*COTB)
        A2P=ZAP*GG1*(TLRM+ZAP*VPVM*GG1/4.0)/2.0
        A2P=A2P +VP2*ZAP*GG1**2/4.0 +4.0*FT2*EP2
        A3=SR2*ZAP*GG1**2/8.0
        A4=-VVP*ZAP*GG1**2/4.0 +4.0*FT2*EP*ATR
C        
        MT1R=SR2*FT*MT-SR2*V*ZAP*GGP/8.0 +A1/(4.0*A0)
        MT2R=SR2*FT*MT-SR2*V*ZAP*GGP/8.0 -A1/(4.0*A0)
        MT1P=SR2*VP*ZAP*GGP/8.0 +A1P/(4.0*A0)
        MT2P=SR2*VP*ZAP*GGP/8.0 -A1P/(4.0*A0)
        MT1RR=FT2 -ZAP*GGP/8.0 -A1**2/(8.0*A0**3) +A2/(4.0*A0)
        MT2RR=FT2 -ZAP*GGP/8.0 +A1**2/(8.0*A0**3) -A2/(4.0*A0)
        MT1PP=ZAP*GGP/8.0 -A1P**2/(8.0*A0**3) +A2P/(4.0*A0)
        MT2PP=ZAP*GGP/8.0 +A1P**2/(8.0*A0**3) -A2P/(4.0*A0)
        MT1RRR=3.0*A1**3/(16.0*A0**3) 
        MT1RRR=MT1RRR/(A0**2) -3.0*A1*A2/(8.0*A0**3) 
     $                            +3.0*V*A3/(4.0*A0)
        MT2RRR=-MT1RRR
        MT1PPP=3.0*A1P**3/(16.0*A0**3) 
        MT1PPP=MT1PPP/(A0**2) -3.0*A1P*A2P/(8.0*A0**3) 
     $                             +3.0*VP*A3/(4.0*A0)
        MT2PPP=-MT1PPP
        MT1RP=-A1*A1P/(8.0*A0**3) +A4/(4.0*A0)
        MT2RP=-MT1RP
        MT1PRR=3.0*A1P*A1**2/(16.0*A0**3)
        MT1PRR=MT1PRR/(A0**2)
     $          -(A2*A1P+2.0*A1*A4)/(8.0*A0**3) -VP*A3/(4.0*A0)
        MT2PRR=-MT1PRR
        MT1RPP=3.0*A1*A1P**2/(16.0*A0**3)
        MT1RPP=MT1RPP/(A0**2)
     $          -(A1*A2P+2.0*A1P*A4)/(8.0*A0**3) -V*A3/(4.0*A0)
        MT2RPP=-MT1RPP
      ELSEIF(RTT.EQ.0.0) THEN
        IF(INRAD.EQ.2.OR.TANB.EQ.1.0) THEN
          IF(EP.EQ.0.0.AND.TLRM.EQ.0.0) THEN
            IF(ATR.EQ.0.0.AND.ATI.EQ.0.0) THEN
              ISPECT=1
              MT1R=SR2*V*FT2
              MT2R=SR2*V*FT2
              MT1P=0.0
              MT2P=0.0
              MT1RR=FT2 
              MT2RR=FT2 
              MT1PP=0.0
              MT2PP=0.0
              MT1RRR=0.0
              MT2RRR=0.0
              MT1PPP=0.0
              MT2PPP=0.0
              MT1RP=0.0
              MT2RP=0.0
              MT1PRR=0.0
              MT2PRR=0.0
              MT1RPP=0.0
              MT2RPP=0.0
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(RTT.NE.0.0 .OR. ISPECT.EQ.1) THEN
        SQVT1=2.0*(3.0*MT1R*MT1RR+MST1SQ*MT1RRR)
        SQVT1=SQVT1*LOG(MST1SQ/QQQ2)
        SQVT1=SQVT1 +2.0*MT1R**3/MST1SQ +9.0*MT1R*MT1RR
        SQVT1=SQVT1+MST1SQ*MT1RRR
        SQVT2=2.0*(3.0*MT2R*MT2RR+MST2SQ*MT2RRR)
        SQVT2=SQVT2*LOG(MST2SQ/QQQ2)
        SQVT2=SQVT2 +2.0*MT2R**3/MST2SQ +9.0*MT2R*MT2RR
        SQVT2=SQVT2+MST2SQ*MT2RRR
        SQVRRR=SQVT1+SQVT2
C          
        SQVT1=2.0*(3.0*MT1P*MT1PP+MST1SQ*MT1PPP)
        SQVT1=SQVT1*LOG(MST1SQ/QQQ2)
        SQVT1=SQVT1 +2.0*MT1P**3/MST1SQ + 9.0*MT1P*MT1PP
        SQVT1=SQVT1+MST1SQ*MT1PPP
        SQVT2=2.0*(3.0*MT2P*MT2PP+MST2SQ*MT2PPP)
        SQVT2=SQVT2*LOG(MST2SQ/QQQ2)
        SQVT2=SQVT2 +2.0*MT2P**3/MST2SQ +9.0*MT2P*MT2PP
        SQVT2=SQVT2 +MST2SQ*MT2PPP
        SQVPPP = SQVT1 + SQVT2
C
        SQVT1=2.0*MT1R*MT1RP+MT1P*MT1RR+MST1SQ*MT1PRR
        SQVT1=2.0*SQVT1*LOG(MST1SQ/QQQ2)
        SQVT1=SQVT1 +2.0*MT1P*MT1R**2/MST1SQ 
        SQVT1=SQVT1+3.0*MT1P*MT1RR+6.0*MT1R*MT1RP
        SQVT1=SQVT1+MST1SQ*MT1PRR
        SQVT2=2.0*MT2R*MT2RP+MT2P*MT2RR+MST2SQ*MT2PRR
        SQVT2=2.0*SQVT2*LOG(MST2SQ/QQQ2)
        SQVT2=SQVT2 +2.0*MT2P*MT2R**2/MST2SQ 
        SQVT2=SQVT2+3.0*MT2P*MT2RR+6.0*MT2R*MT2RP
        SQVT2=SQVT2+MST2SQ*MT2PRR
        SQVPRR=SQVT1+SQVT2
C
        SQVT1=2.0*MT1P*MT1RP+MT1R*MT1PP+MST1SQ*MT1RPP
        SQVT1=2.0*SQVT1*LOG(MST1SQ/QQQ2)
        SQVT1=SQVT1 +2.0*MT1R*MT1P**2/MST1SQ 
        SQVT1=SQVT1+3.0*MT1R*MT1PP+6.0*MT1P*MT1RP
        SQVT1=SQVT1+MST1SQ*MT1RPP
        SQVT2=2.0*MT2P*MT2RP+MT2R*MT2PP+MST2SQ*MT2RPP
        SQVT2=2.0*SQVT2*LOG(MST2SQ/QQQ2)
        SQVT2=SQVT2 +2.0*MT2R*MT2P**2/MST2SQ 
        SQVT2=SQVT2+3.0*MT2R*MT2PP+6.0*MT2P*MT2RP
        SQVT2=SQVT2+MST2SQ*MT2RPP
        SQVRPP=SQVT1+SQVT2
C
        FVRRR=-2.0*SR2*FT4*V*(6.0*LOG(MT2/QQQ2) + 13.0)
      ENDIF
C
      IF(RTT.EQ.0.0 .AND. ISPECT.EQ.0) THEN
        ALPHAT=(MSTL2 + MSTR2)/2.0 + MT2
        ALPHAT=ALPHAT +VP2*(1.0-TANB**2)*ZAP*GGP/8.0
        GGP1SQ= ZAP*GGP**2 +ZAP*GG1**2
C
        SQVRRR=12.0*FT4*LOG(ALPHAT/MT2)
        TEMPSQ=-FT2*ZAP*GGP +GGP1SQ/16.0
        SQVRRR=SQVRRR +3.0*TEMPSQ*LOG(ALPHAT/QQQ2)
        SQVRRR=SQVRRR -8.0*FT4 -9.0*FT2*ZAP*GGP/2.0
        SQVRRR=SQVRRR +9.0*GGP1SQ/32.0
        TEMPSQ=8.0*V2*(FT2-ZAP*GGP/8.0)**2
        TEMPSQ=TEMPSQ +3.0*V2*ZAP*GG1**2/8.0
        TEMPSQ=TEMPSQ +6.0*FT2*EP2*COTB**2
        SQVRRR=SQVRRR +TEMPSQ*(FT2-ZAP*GGP/8.0)/ALPHAT
        SQVRRR=SQVRRR*SR2*V
C
        SQVPPP=3.0*GGP1SQ*(2.0*LOG(ALPHAT/QQQ2)+3.0)/32.0
        TEMPSQ=ZAP*GGP*(ZAP*GGP**2+3.0*GG1**2)*VP2/ALPHAT/64.0
        SQVPPP=SQVPPP+TEMPSQ
        TEMPSQ=3.0*FT2*EP2*ZAP*GGP/ALPHAT/4.0
        SQVPPP=(SQVPPP+TEMPSQ)*SR2*VP
C
        TEMPSQ=FT2*ZAP*GGP -GGP1SQ/8.0
        SQVPRR=TEMPSQ*(2.0*LOG(ALPHAT/QQQ2)+3.0)
        TEMPSQ=4.0*ZAP*GGP*(FT2-ZAP*GGP/8.0)-FT2*ZAP*GG1**2
        TEMPSQ=TEMPSQ +3.0*ZAP*GGP*GG1**2/16.0
        TEMPSQ=V2*TEMPSQ+EP2*FT2*ZAP*GGP*(2.0+COTB**2)
        TEMPSQ=TEMPSQ-16.0*EP2*FT4
        SQVPRR=(SQVPRR+TEMPSQ/ALPHAT)*SR2*VP/4.0
C
        TEMPSQ=FT2*ZAP*GGP -GGP1SQ/8.0
        SQVRPP=TEMPSQ*(2.0*LOG(ALPHAT/QQQ2)+3.0)
        TEMPSQ=GGP1SQ*(FT2-ZAP*GGP/8.0)-ZAP*GGP*GG1**2/4.0
        TEMPSQ=VP2*TEMPSQ/2.0 +8.0*EP2*FT4
        TEMPSQ=TEMPSQ+EP2*FT2*ZAP*GGP*(1.0+2.0*COTB**2)
        SQVRPP=(SQVRPP+TEMPSQ/ALPHAT)*SR2*V/4.0
C
        FVRRR = 0.0
C      
C        Fermion part (FRRR) is already combined 
C        with the squark part. 
C
      ENDIF     
C
C
C                  SBOTTOM TERMS
C
      ISPECB=0          
      RBB=(BLRM-VPVM*ZAP*GG2/4.0)**2
     $      +4.0*MB2*(EP*TANB+ABR)**2+4.0*MB2*ABI**2
C
      IF(RBB.GT.0.0) THEN      
        B0=SQRT(RBB)
        B1=V*ZAP*GG2*(BLRM-ZAP*VPVM*GG2/4.0)/SR2
        B1=B1+4.0*SR2*FB*MB*EP*(ABR+EP*TANB)
        B2=ZAP*GG2*(BLRM-ZAP*VPVM*GG2/4.0)/2.0
        B2=B2 +V2*ZAP*GG2**2/4.0 +4.0*FB2*EP2
        B1P=-VP*ZAP*GG2*(BLRM-ZAP*VPVM*GG2/4.0)/SR2
        B1P=B1P+4.0*SR2*FB*MB*(AB2+EP*ABR*TANB)
        B2P=-ZAP*GG2*(BLRM-ZAP*VPVM*GG2/4.0)/2.0
        B2P=B2P +VP2*ZAP*GG2**2/4.0 +4.0*FB2*AB2
        B3=SR2*ZAP*GG2**2/8.0
        B4=-VVP*ZAP*GG2**2/4.0 +4.0*FB2*EP*ABR
C
        MB1R=SR2*V*ZAP*GGP/8.0 +B1/(4.0*B0)
        MB2R=SR2*V*ZAP*GGP/8.0 -B1/(4.0*B0)
        MB1P=SR2*FB*MB -SR2*VP*ZAP*GGP/8.0 +B1P/(4.0*B0)
        MB2P=SR2*FB*MB -SR2*VP*ZAP*GGP/8.0 -B1P/(4.0*B0)
        MB1RR=ZAP*GGP/8.0 -B1**2/(8.0*B0**3) +B2/(4.0*B0)
        MB2RR=ZAP*GGP/8.0 +B1**2/(8.0*B0**3) -B2/(4.0*B0)
        MB1PP=FB2 -ZAP*GGP/8.0 
        MB2PP=MB1PP 
        MB1PP=MB1PP -B1P**2/(8.0*B0**3) +B2P/(4.0*B0)
        MB2PP=MB2PP +B1P**2/(8.0*B0**3) -B2P/(4.0*B0)
        MB1RRR=3.0*B1**3/(16.0*B0**3) 
        MB1RRR=MB1RRR/(B0**2) -3.0*B1*B2/(8.0*B0**3) 
     $                            +3.0*V*B3/(4.0*B0)
        MB2RRR=-MB1RRR
        MB1PPP=3.0*B1P**3/(16.0*B0**3) 
        MB1PPP=MB1PPP/(B0**2) -3.0*B1P*B2P/(8.0*B0**3) 
        MB1PPP=MB1PPP +3.0*VP*B3/(4.0*B0)
        MB2PPP=-MB1PPP 
        MB1RP=-B1*B1P/(8.0*B0**3) +B4/(4.0*B0)
        MB2RP=-MB1RP
        MB1PRR=3.0*B1P*B1**2/(16.0*B0**3)
        MB1PRR=MB1PRR/(B0**2) -(B2*B1P+2.0*B1*B4)/(8.0*B0**3)
        MB1PRR=MB1PRR -VP*B3/(4.0*B0)
        MB2PRR=-MB1PRR
        MB1RPP=3.0*B1*B1P**2/(16.0*B0**3)
        MB1RPP=MB1RPP/(B0**2) -(B1*B2P+2.0*B1P*B4)/(8.0*B0**3)
        MB1RPP=MB1RPP -V*B3/(4.0*B0)
        MB2RPP=-MB1RPP
      ELSEIF(RBB.EQ.0.0) THEN
        IF(INRAD.EQ.2.OR.TANB.EQ.1.0) THEN
          IF(EP.EQ.0.0.AND.BLRM.EQ.0.0) THEN
            IF(ABR.EQ.0.0.AND.ABI.EQ.0.0) THEN
              ISPECB=1
              MB1R=0.0
              MB2R=0.0
              MB1P=SR2*VP*FB2 
              MB2P=SR2*VP*FB2
              MB1RR=0.0
              MB2RR=0.0
              MB1PP=FB2
              MB2PP=FB2
              MB1RRR=0.0
              MB2RRR=0.0
              MB1PPP=0.0
              MB2PPP=0.0
              MB1RP=0.0
              MB2RP=0.0
              MB1PRR=0.0
              MB1PRR=0.0
              MB2PRR=0.0
              MB1RPP=0.0
              MB2RPP=0.0
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C
      IF(RBB.NE.0.0 .OR. ISPECB.EQ.1) THEN
        SQVB1=2.0*(3.0*MB1R*MB1RR+MSB1SQ*MB1RRR)
        SQVB1=SQVB1*LOG(MSB1SQ/QQQ2)
        SQVB1=SQVB1 +2.0*MB1R**3/MSB1SQ +9.0*MB1R*MB1RR
        SQVB1=SQVB1+MSB1SQ*MB1RRR
        SQVB2=2.0*(3.0*MB2R*MB2RR+MSB2SQ*MB2RRR)
        SQVB2=SQVB2*LOG(MSB2SQ/QQQ2)
        SQVB2=SQVB2 +2.0*MB2R**3/MSB2SQ +9.0*MB2R*MB2RR
        SQVB2=SQVB2+MSB2SQ*MB2RRR
        SQVRRR = SQVRRR + SQVB1 + SQVB2
C
        SQVB1=2.0*(3.0*MB1P*MB1PP+MSB1SQ*MB1PPP)
        SQVB1=SQVB1*LOG(MSB1SQ/QQQ2)
        SQVB1=SQVB1 +2.0*MB1P**3/MSB1SQ +9.0*MB1P*MB1PP
        SQVB1=SQVB1+MSB1SQ*MB1PPP
        SQVB2=2.0*(3.0*MB2P*MB2PP+MSB2SQ*MB2PPP)
        SQVB2=SQVB2*LOG(MSB2SQ/QQQ2)
        SQVB2=SQVB2 +2.0*MB2P**3/MSB2SQ +9.0*MB2P*MB2PP
        SQVB2=SQVB2+MSB2SQ*MB2PPP
        SQVPPP= SQVPPP+SQVB1+SQVB2
C
        SQVB1=2.0*MB1R*MB1RP+MB1P*MB1RR+MSB1SQ*MB1PRR
        SQVB1=2.0*SQVB1*LOG(MSB1SQ/QQQ2)
        SQVB1=SQVB1 +2.0*MB1P*MB1R**2/MSB1SQ 
        SQVB1=SQVB1 +3.0*MB1P*MB1RR +6.0*MB1R*MB1RP
        SQVB1=SQVB1+MSB1SQ*MB1PRR
        SQVB2=2.0*MB2R*MB2RP+MB2P*MB2RR+MSB2SQ*MB2PRR
        SQVB2=2.0*SQVB2*LOG(MSB2SQ/QQQ2)
        SQVB2=SQVB2 +2.0*MB2P*MB2R**2/MSB2SQ 
        SQVB2=SQVB2 +3.0*MB2P*MB2RR +6.0*MB2R*MB2RP
        SQVB2=SQVB2+MSB2SQ*MB2PRR
        SQVPRR=SQVPRR+SQVB1+SQVB2
C
        SQVB1=2.0*MB1P*MB1RP+MB1R*MB1PP+MSB1SQ*MB1RPP
        SQVB1=2.0*SQVB1*LOG(MSB1SQ/QQQ2)
        SQVB1=SQVB1 +2.0*MB1R*MB1P**2/MSB1SQ 
        SQVB1=SQVB1+3.0*MB1R*MB1PP+6.0*MB1P*MB1RP
        SQVB1=SQVB1+MSB1SQ*MB1RPP
        SQVB2=2.0*MB2P*MB2RP+MB2R*MB2PP+MSB2SQ*MB2RPP
        SQVB2=2.0*SQVB2*LOG(MSB2SQ/QQQ2)
        SQVB2=SQVB2 +2.0*MB2R*MB2P**2/MSB2SQ 
        SQVB2=SQVB2 +3.0*MB2R*MB2PP +6.0*MB2P*MB2RP
        SQVB2=SQVB2+MSB2SQ*MB2RPP
        SQVRPP=SQVRPP+SQVB1+SQVB2
C
        IF(MB2.EQ.0.0) THEN
          FVPPP=0.0
        ELSE IF(MB2.NE.0.0) THEN
          FVPPP=-2.0*SR2*FB4*VP*(6.0*LOG(MB2/QQQ2)+13.0)
        ENDIF
C
      ENDIF
C
      IF(RBB.EQ.0.0 .AND. ISPECB.EQ.0) THEN
        ALPHAB=(MSBL2+MSBR2)/2.0 +MB2 
        ALPHAB=ALPHAB -VP2*(1.0-TANB**2)*ZAP*GGP/8.0
        GGP2SQ=ZAP*GGP**2 +ZAP*GG2**2
C
        BSQ=3.0*GGP2SQ*(2.0*LOG(ALPHAB/QQQ2)+3.0)/8.0
        TEMPSQ=V2*(ZAP*GGP**2 +3.0*ZAP*GG2**2)/16.0 
     $                                   +3.0*FB2*EP2
        BSQ=(BSQ +ZAP*GGP*TEMPSQ/ALPHAB)*SR2*V/4.0
        SQVRRR=SQVRRR+BSQ
C
        BSQ=12.0*FB4*LOG(ALPHAB/MB2) -8.0*FB4
        TEMPSQ=-FB2*ZAP*GGP +GGP2SQ/16.0
        BSQ=BSQ+3.0*TEMPSQ*(LOG(ALPHAB/QQQ2)+1.5)
        TEMPSQ=8.0*VP2*(FB2-ZAP*GGP/8.0)**2 
     $     +3.0*VP2*ZAP*GG2**2/8.0 +6.0*FB2*EP2*TANB**2
        BSQ=BSQ +(FB2-ZAP*GGP/8.0)*TEMPSQ/ALPHAB
        BSQ=BSQ*SR2*VP
        SQVPPP=SQVPPP+BSQ
C
        TEMPSQ=0.5*(FB2*ZAP*GGP -GGP2SQ/8.0)
        BSQ=TEMPSQ*(LOG(ALPHAB/QQQ2)+1.5)
        TEMPSQ=(FB2 -ZAP*GGP/8.0)*GGP2SQ -ZAP*GGP*GG2**2/4.0
        TEMPSQ=V2*TEMPSQ/4.0 +4.0*FB4*EP2 -FB2*EP2*ZAP*GGP/2.0
        TEMPSQ=(TEMPSQ-FB2*EP2*ZAP*GGP*TANB**2)/ALPHAB/2.0
        BSQ=(BSQ+TEMPSQ)*SR2*VP
        SQVPRR=SQVPRR+BSQ
C
        TEMPSQ=0.5*(FB2*ZAP*GGP -GGP2SQ/8.0)
        BSQ=TEMPSQ*(LOG(ALPHAB/QQQ2)+1.5)
        TEMPSQ=4.0*ZAP*GGP*(FB2 -ZAP*GGP/8.0)**2 
     $          -FB2*ZAP*GG2**2 +3.0*ZAP*GGP*GG2**2/16.0
        TEMPSQ=VP2*TEMPSQ-16.0*FB4*EP2 
        TEMPSQ=TEMPSQ+FB2*EP2*ZAP*GGP*(TANB**2 +0.5)
        BSQ=(BSQ +TEMPSQ/ALPHAB/4.0)*SR2*V
        SQVRPP=SQVRPP+BSQ
C
        FVPPP=0.0
C
C        Fermion part (FPPP) is already combined 
C        with the squark part. 
C
      ENDIF
C
C
      VRRR=3.0*(SQVRRR+FVRRR)/(32.0*PI2)
      VRRR=VRRR/6.0
C
      VPPP=3.0*(SQVPPP+FVPPP)/(32.0*PI2)
      VPPP=VPPP/6.0
C
      VPRR=3.0*(SQVPRR)/(32.0*PI2)
      VPRR=VPRR/2.0
C
      VRPP=3.0*(SQVRPP)/(32.0*PI2)
      VRPP=VRPP/2.0
C
C
C      Note in the following that the angle ALFAH
C      calculated in the subroutine SSMHN must
C      be input.
C
      CA2=(COS(ALFAH))**2
      SA2=(SIN(ALFAH))**2
      DVHLL=-VRRR*CA2*SIN(ALFAH)
      DVHLL=DVHLL +VPRR*(CA2-2.0*SA2)*COS(ALFAH)
      DVHLL=DVHLL +VRPP*(2.0*CA2-SA2)*SIN(ALFAH)
      DVHLL=DVHLL +VPPP*SA2*COS(ALFAH)
C
      DVHLL=3.0*DVHLL
      DVHLL=-DVHLL
C
C      Finally, multiply bt the coefficient of the
C      tree-level Lagrangian level term (COEFF.) 
C      so that the answer may be written as:  
C        DW = (COEFF.)**2 
C            * (TREE-LEVEL ANGULAR DEPENDENCE + DVHLL)
C
C              *(LAMBDA KINEMATIC FCN)**0.5/(8*PI*MHH**3)
C
C
      DVHLL=4.0*SQRT((1.-SN2THW)/G2)*DVHLL/AMZ
C
C
1000  DELHLL=DVHLL
      RETURN
      END

#include "isajet/pilot.h"
C-----------------------------------------------------------------
      SUBROUTINE SUGEFF(G0,SIG1,SIG2)
C-----------------------------------------------------------------
C
C     Compute Higgs mass shift due to 1-loop effective potential
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sugpas.inc"
#include "isajet/sugmg.inc"
      REAL G0(29),SIG1,SIG2
      REAL DT1,DELT1S,SIG1T1,DMSDV2,FB,FT,MST2,MSB2,MSB1,SIG2B1,
     $SIG1B1,SIG2B2,SIG1B2,DB1,SIG1T2,SIG2T1,DELB1S,SIG2T2,MST1,COS2W, 
     $MT,PI,COTB,TANB,MB,SIG2B,SIG1B,G,GGP,SIG2T,E,FAC,SIG1T,QS,
     $BETA,SINB,COSB,FAC4,FL,ML,SIG1L,SIG2L,MSL1,MSL2,
     $DELL1S,DL1,SIG1L1,SIG2L1,SIG1L2,SIG2L2
C
      G=G2
      TANB=XTANB
      COS2W=1.-XW
      COTB=1./TANB
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      PI=4.*ATAN(1.)
      FAC=3./8./PI**2
      FAC4=FAC/3.
      E=EXP(1.)
      QS=HIGFRZ**2
C-----CALCULATE TOP AND BOTTOM CONTRIBUTIONS; USE RUNNING MASSES--
      FL=G0(4)
      FB=G0(5)
      FT=G0(6)
      ML=FL*VEV*COSB
      MT=FT*VEV*SINB
      MB=FB*VEV*COSB
      SIG1T=0.
      SIG2T=-FAC*MT**2*G0(6)**2*LOG(MT**2/E/QS)
      SIG1B=-FAC*MB**2*G0(5)**2*LOG(MB**2/E/QS)
      SIG2B=0.
      SIG1L=-FAC4*ML**2*G0(4)**2*LOG(ML**2/E/QS)
      SIG2L=0.
      GGP=(G**2+GP**2)/2.
      MST1=MSS(12)
      MST2=MSS(13)
      MSB1=MSS(10)
      MSB2=MSS(11)
      MSL1=MSS(21)
      MSL2=MSS(22)
C-----CALCULATE STOP_1 CONTRIBUTION -------------------------------
      DELT1S=(.5*(G0(24)-G0(23))+(8*COS2W-5.)*GGP*
     $        (VP*VP-V*V)/12.)**2+G0(6)**2*V*V*(G0(12)-MU*COTB)**2
      DT1=.5*(G0(24)-G0(23))+(8*COS2W-5.)*GGP*
     $       (VP*VP-V*V)/12.
      DMSDV2=GGP/4.-(2*DT1*(8*COS2W-5.)*GGP/12.-
     $        FT**2*MU*(G0(12)*TANB-MU))/2./SQRT(DELT1S)
      SIG1T1=FAC/2.*MST1**2*LOG(MST1**2/E/QS)*DMSDV2
      DMSDV2=-GGP/4.+FT**2-(-2*DT1*(8*COS2W-5.)*GGP/12.+
     $        FT**2*G0(12)*(G0(12)-MU*COTB))/2./SQRT(DELT1S)
      SIG2T1=FAC/2.*MST1**2*LOG(MST1**2/E/QS)*DMSDV2
C-----CALCULATE STOP_2 CONTRIBUTION -------------------------------
      DMSDV2=GGP/4.+(2*DT1*(8*COS2W-5.)*GGP/12.-
     $        FT**2*MU*(G0(12)*TANB-MU))/2./SQRT(DELT1S)
      SIG1T2=FAC/2.*MST2**2*LOG(MST2**2/E/QS)*DMSDV2
      DMSDV2=-GGP/4.+FT**2+(-2*DT1*(8*COS2W-5.)*GGP/12.+
     $        FT**2*G0(12)*(G0(12)-MU*COTB))/2./SQRT(DELT1S)
      SIG2T2=FAC/2.*MST2**2*LOG(MST2**2/E/QS)*DMSDV2
C-----CALCULATE SBOT_1 CONTRIBUTION -------------------------------
      DELB1S=(.5*(G0(24)-G0(22))-(4*COS2W-1.)*GGP*
     $  (VP*VP-V*V)/12.)**2+G0(5)**2*VP*VP*(G0(11)-MU*TANB)**2
      DB1=.5*(G0(24)-G0(22))-(4*COS2W-1.)*GGP*
     $       (VP*VP-V*V)/12.
      DMSDV2=-GGP/4.+FB**2-(-2*DB1*(4*COS2W-1.)*GGP/12.+
     $        FB**2*G0(11)*(G0(11)-MU*TANB))/2./SQRT(DELB1S)
      SIG1B1=FAC/2.*MSB1**2*LOG(MSB1**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.-(2*DB1*(4*COS2W-1.)*GGP/12.-
     $        FB**2*MU*(G0(11)*COTB-MU))/2./SQRT(DELB1S)
      SIG2B1=FAC/2.*MSB1**2*LOG(MSB1**2/E/QS)*DMSDV2
C-----CALCULATE SBOT_2 CONTRIBUTION -------------------------------
      DMSDV2=-GGP/4.+FB**2+(-2*DB1*(4*COS2W-1.)*GGP/12.+
     $        FB**2*G0(11)*(G0(11)-MU*TANB))/2./SQRT(DELB1S)
      SIG1B2=FAC/2.*MSB2**2*LOG(MSB2**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.+(2*DB1*(4*COS2W-1.)*GGP/12.-
     $        FB**2*MU*(G0(11)*COTB-MU))/2./SQRT(DELB1S)
      SIG2B2=FAC/2.*MSB2**2*LOG(MSB2**2/E/QS)*DMSDV2
C-----CALCULATE STAU_1 CONTRIBUTION -------------------------------
      DELL1S=(.5*(G0(21)-G0(20))-(4*COS2W-3.)*GGP*
     $  (VP*VP-V*V)/4.)**2+G0(4)**2*VP*VP*(G0(10)-MU*TANB)**2
      DL1=.5*(G0(21)-G0(20))-(4*COS2W-3.)*GGP*
     $       (VP*VP-V*V)/4.
      DMSDV2=-GGP/4.+FL**2-(-2*DL1*(4*COS2W-3.)*GGP/4.+
     $        FL**2*G0(10)*(G0(10)-MU*TANB))/2./SQRT(DELL1S)
      SIG1L1=FAC4/2.*MSL1**2*LOG(MSL1**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.-(2*DL1*(4*COS2W-3.)*GGP/4.-
     $        FL**2*MU*(G0(10)*COTB-MU))/2./SQRT(DELL1S)
      SIG2L1=FAC4/2.*MSL1**2*LOG(MSL1**2/E/QS)*DMSDV2
C-----CALCULATE STAU_2 CONTRIBUTION -------------------------------
      DMSDV2=-GGP/4.+FL**2+(-2*DL1*(4*COS2W-3.)*GGP/4.+
     $        FL**2*G0(10)*(G0(10)-MU*TANB))/2./SQRT(DELL1S)
      SIG1L2=FAC4/2.*MSL2**2*LOG(MSL2**2/E/QS)*DMSDV2
      DMSDV2=GGP/4.+(2*DL1*(4*COS2W-3.)*GGP/4.-
     $        FL**2*MU*(G0(10)*COTB-MU))/2./SQRT(DELL1S)
      SIG2L2=FAC4/2.*MSL2**2*LOG(MSL2**2/E/QS)*DMSDV2
C-----ADD ALL TERMS ------------------------------------------------
      SIG1=SIG1B+SIG1B1+SIG1B2+SIG1T+SIG1T1+SIG1T2+
     $SIG1L+SIG1L1+SIG1L2
      SIG2=SIG2B+SIG2B1+SIG2B2+SIG2T+SIG2T1+SIG2T2+
     $SIG2L+SIG2L1+SIG2L2
      RETURN
      END

#include "isajet/pilot.h"
      SUBROUTINE SSHSF
C-----------------------------------------------------------------------
C
C          Calculates the partial decay widths of 
C          the Higgs bosons into sfermions.
C          calculated by X. Tata
C          program by M. Bisset
C
C          10/23/93: modified by H. Baer, 10/8/96
C          Intra-flavor sfermion mixing is neglected
C          for all flavors EXCEPT for stops, sbottoms and staus.
C
C
C         10/23/93
C          It is assumed that the A-terms are real.
C          In addition, all coefficients of the sfermion
C          trilinear terms from the superpotential 
C          EXCEPT the stop (AAT), sbottom (AAB) and stau (AAL)
C          coefficients are set to zero.
C
C     ===> Code for the general case removing all these 
C          artificial restrictions is present below.  
C          The preceeding restrictions are specified
C          by giving special values to some variables
C          This is discussed in two sections beginning
C          with the symbols  (*@&*)  in the code below.
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/sstype.inc"
C
C
      REAL SR2,PI,GG,TW2,BETA,DSA,DCA,DSB,DCB,MH
      REAL EP,TANB,COTB,ATERM,MSFMIX,THETSF,SIN2B
      REAL TEMP,TEMP1,TEMP2,YA1,YA2
      REAL SINA,COSA,SINA2,COSA2,M1,M2,M12,LAMB
      REAL SINAU,COSAU,SINAD,COSAD
      REAL A11,A22,A12,B11,B22,B12,C11,C12,C21,C22
      REAL ASQ,BSQ,CSQ,DWSF
      REAL DWSFL,DWSFH,DWSFP,DWSFC,SSXLAM
      REAL ASMB,MBMB,MBQ,ASMT,MTMT,MTQ,SUALFS
      DOUBLE PRECISION SSMQCD
      DIMENSION ATERM(12),MSFMIX(12,2),THETSF(12)
      DIMENSION ASQ(10,3),BSQ(9),CSQ(6,4)
      DIMENSION DWSF(12,4),DWSFL(12,4),DWSFH(12,4)
      DIMENSION DWSFP(12,4),DWSFC(6,4)
      INTEGER II,IJ,JJ,IC,IJU,IJD,NUMH
C
C
      SR2=SQRT(2.0)
      PI=4.0*ATAN(1.0)
      TW2=SN2THW/(1.0-SN2THW)
      GG=SQRT(4.0*PI*ALFAEM/SN2THW)
      EP=TWOM1        
C
      TANB=1.0/RV2V1
      COTB=RV2V1
      BETA=ATAN(1.0/RV2V1)
      DSA=SIN(ALFAH)
      DCA=COS(ALFAH)
      DSB=SIN(BETA)
      DCB=COS(BETA)
      SIN2B=2.0*DSB*DCB
C      
C      Set A-terms.
C      (all A-terms are assumed to be real)
C      The A-terms are loaded into the array ATERM(12)
C      in the following way: 
C          ATERM(1)=selectron A-term
C          ATERM(2)=smuon A-term
C          ATERM(3)=stau A-term
C          ATERM(4)=up squark A-term
C          ATERM(5)=charm squark A-term
C          ATERM(6)=down squark A-term
C          ATERM(7)=strange squark A-term
C          ATERM(8)=sbottom A-term
C          ATERM(9)=stop A-term
C          ATERM(10)=selectronic sneutrino A-term
C          ATERM(11)=smuonic sneutrino A-term
C          ATERM(12)=stauonic sneutrino A-term
C
      DO 10 II=1,7
        ATERM(II)=0.0
10    CONTINUE
      ATERM(3)=AAL
      ATERM(8)=AAB
      ATERM(9)=AAT  
      DO 20 II=10,12
        ATERM(II)=0.0
20    CONTINUE
C
C      Set mixing parameters.
C      The intra-flavor-mixed sfermion masses are loaded into
C      the array MSFMIX(12,2) where (#,1) is the lighter 
C      mixed sfermion mass of a given flavor and (#,2) is the 
C      heavier sfermion mass.  The sfermionic mixing angles are
C      loaded into the array THETSF(12).  The identities of the
C      elements of these arrays are given below:
C        MSFMIX(1,*)=mixed selectron masses  
C                                 THETSF(1)=selectron mixing angle
C        MSFMIX(2,*)=mixed smuon masses   
C                                 THETSF(2)=smuon mixing angle
C        MSFMIX(3,*)=mixed stau masses  
C                                 THETSF(3)=stau mixing angle
C        MSFMIX(4,*)=mixed up squark masses 
C                                 THETSF(4)=up squark mixing angle
C        MSFMIX(5,*)=mixed charm squark masses 
C                                 THETSF(5)=charm squark mixing angle
C        MSFMIX(6,*)=mixed down squark masses 
C                                 THETSF(6)=down squark mixing angle
C        MSFMIX(7,*)=mixed strange squark masses  
C                                 THETSF(7)=strange squark mixing angle
C        MSFMIX(8,*)=mixed sbottom masses 
C                                 THETSF(8)=sbottom mixing angle
C        MSFMIX(9,*)=mixed stop masses  
C                                 THETSF(9)=stop mixing angle
C        For sneuterinos MSFMIX(#,2)=0.0, THETSF(#)=0.0 ; #=10-12
C        Yukawa contributions from D-terms to the sneutrino masses
C        are supposed to be added in here. 
C        MSFMIX(10,1)= selectronic sneutrino mass with D-terms 
C        MSFMIX(11,1)= smuonic sneutrino mass with D-terms
C        MSFMIX(12,1)= stauonic sneutrino mass with D-terms
C
      DO 30 II=10,12
        MSFMIX(II,2)=0.0
        THETSF(II)=0.0
30    CONTINUE
C
C
C       (*@&*) 10/24/93 - Special conditions used ---
C        set all mixing angles EXCEPT stop, sbottom, stau to zero.
C        For all EXCEPT st, sb and stau, set mixed sfermion masses
C        to bare sfermion masses:
C            MSFMIX(#,1) = Left sfermion mass
C            MSFMIX(#,2) = Right sfermion mass  ;  # = 1-8
C         but
C            MSFMIX(9,1) = AMT1SS  
C            MSFMIX(9,2) = AMT2SS              , etc.
C
C        (The choice of which to call Left and which to call
C         Right is based on the definition of the sfermion
C         mixing angle theta_sf :
C            sfermion_1 = cos(theta_sf) * sfermion_L
C                            - sin(theta_sf) * sfermion_R
C            sfermion_2 = sin(theta_sf) * sfermion_L
C                            + cos(theta_sf) * sfermion_R
C          Thus if we set theta_sf = 0, then
C                 sfermion_1 = sfermion_L
C           and   sfermion_2 = sfermion_R .               )
C
      DO 40 II=1,7
        THETSF(II)=0.0
40    CONTINUE
      MSFMIX(1,1)=AMELSS
      MSFMIX(1,2)=AMERSS
      MSFMIX(2,1)=AMMLSS
      MSFMIX(2,2)=AMMRSS
      MSFMIX(3,1)=AML1SS
      MSFMIX(3,2)=AML2SS
      THETSF(3)=THETAL
      MSFMIX(4,1)=AMULSS
      MSFMIX(4,2)=AMURSS
      MSFMIX(5,1)=AMCLSS
      MSFMIX(5,2)=AMCRSS
      MSFMIX(6,1)=AMDLSS
      MSFMIX(6,2)=AMDRSS
      MSFMIX(7,1)=AMSLSS
      MSFMIX(7,2)=AMSRSS
      MSFMIX(8,1)=AMB1SS
      MSFMIX(8,2)=AMB2SS
      THETSF(8)=THETAB
      MSFMIX(9,1)=AMT1SS
      MSFMIX(9,2)=AMT2SS
      THETSF(9)=THETAT
      MSFMIX(10,1)=AMN1SS
      MSFMIX(11,1)=AMN2SS
      MSFMIX(12,1)=AMN3SS
C
      DO 1000 NUMH=1,4
        IF(NUMH.EQ.1) THEN
          MH=AMHL
        ELSE IF(NUMH.EQ.2) THEN
          MH=AMHH
        ELSE IF(NUMH.EQ.3) THEN
          MH=AMHA
          GO TO 233
        ELSE IF(NUMH.EQ.4) THEN
          MH=AMHC
          GO TO 333
        ENDIF
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(MH))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(MH))
 
C
C         Scalar neutral Higgses --> sfermions 
C          partial decay widths
C
        IF(NUMH.EQ.1) THEN
          TEMP=GG*AMW*SIN(BETA-ALFAH)/2.0
          YA1=DCA
          YA2=DSA
        ELSE IF(NUMH.EQ.2) THEN
          TEMP=-GG*AMW*COS(BETA-ALFAH)/2.0
          YA1=-DSA
          YA2=DCA
        ENDIF
C
        TEMP1=TEMP*(1.0-TW2/3.0)
        TEMP2=GG*YA1/(AMW*DSB)
        ASQ(4,1)=TEMP1-TEMP2*AMUP**2
        ASQ(5,1)=TEMP1-TEMP2*AMCH**2
        ASQ(9,1)=TEMP1-TEMP2*MTQ**2
C
        TEMP1=-TEMP*(1.0+TW2/3.0)
        TEMP2=GG*YA2/(AMW*DCB)
        ASQ(6,1)=-TEMP1-TEMP2*AMDN**2
        ASQ(7,1)=-TEMP1-TEMP2*AMST**2
        ASQ(8,1)=-TEMP1-TEMP2*MBQ**2
C
        ASQ(10,1)=TEMP*(1.0+TW2)
        TEMP1=TEMP*(TW2-1.0)
        TEMP2=GG*YA2/(AMW*DCB)
        ASQ(1,1)=TEMP1-TEMP2*AME**2
        ASQ(2,1)=TEMP1-TEMP2*AMMU**2
        ASQ(3,1)=TEMP1-TEMP2*AMTAU**2
C
        TEMP1=4.0*TEMP*TW2/3.0
        TEMP2=GG*YA1/(AMW*DSB)
        ASQ(4,2)=TEMP1-TEMP2*AMUP**2
        ASQ(5,2)=TEMP1-TEMP2*AMCH**2
        ASQ(9,2)=TEMP1-TEMP2*MTQ**2
C
        TEMP1=-2.0*TEMP*TW2/3.0
        TEMP2=GG*YA2/(AMW*DCB)
        ASQ(6,2)=TEMP1-TEMP2*AMDN**2
        ASQ(7,2)=TEMP1-TEMP2*AMST**2
        ASQ(8,2)=TEMP1-TEMP2*MBQ**2
C
        ASQ(10,2)=0.0
        TEMP1=-2.0*TEMP*TW2
        TEMP2=GG*YA2/(AMW*DCB)
        ASQ(1,2)=TEMP1-TEMP2*AME**2
        ASQ(2,2)=TEMP1-TEMP2*AMMU**2
        ASQ(3,2)=TEMP1-TEMP2*AMTAU**2
C
        TEMP1=GG/(2.0*AMW*DSB)
        ASQ(4,3)=(EP*YA2 + ATERM(4)*YA1)*TEMP1*AMUP
        ASQ(5,3)=(EP*YA2 + ATERM(5)*YA1)*TEMP1*AMCH 
        ASQ(9,3)=(EP*YA2 + ATERM(9)*YA1)*TEMP1*MTQ 
C
        TEMP1=GG/(2.0*AMW*DCB)
        ASQ(6,3)=(ATERM(6)*YA2 + EP*YA1)*TEMP1*AMDN
        ASQ(7,3)=(ATERM(7)*YA2 + EP*YA1)*TEMP1*AMST
        ASQ(8,3)=(ATERM(8)*YA2 + EP*YA1)*TEMP1*MBQ
C
        ASQ(10,3)=0.0
        ASQ(1,3)=(ATERM(1)*YA2 + EP*YA1)*TEMP1*AME
        ASQ(2,3)=(ATERM(2)*YA2 + EP*YA1)*TEMP1*AMMU
        ASQ(3,3)=(ATERM(3)*YA2 + EP*YA1)*TEMP1*AMTAU
C
C
        DO 150 IJ=1,9
          IF(IJ.LT.4) THEN
            TEMP1=1.0/(16.0*PI*MH**3)
          ELSE 
            TEMP1=3.0/(16.0*PI*MH**3)
          ENDIF
          SINA=SIN(THETSF(IJ))
          COSA=COS(THETSF(IJ))
          SINA2=SINA**2
          COSA2=COSA**2
          M1=MSFMIX(IJ,1)
          M2=MSFMIX(IJ,1)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            A11=ASQ(IJ,1)*COSA2+ASQ(IJ,2)*SINA2 
     $              -2.0*ASQ(IJ,3)*SINA*COSA
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSF(IJ,1)=TEMP1*SQRT(LAMB)*A11**2
          ELSE IF(MH.LE.M12) THEN
            DWSF(IJ,1)=0.0
          ENDIF
C
          M1=MSFMIX(IJ,2)
          M2=MSFMIX(IJ,2)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            A22=ASQ(IJ,1)*SINA2+ASQ(IJ,2)*COSA2 
     $                 +2.0*ASQ(IJ,3)*SINA*COSA
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSF(IJ,2)=TEMP1*SQRT(LAMB)*A22**2
          ELSE IF(MH.LE.M12) THEN
            DWSF(IJ,2)=0.0
          ENDIF
C          
          M1=MSFMIX(IJ,1)
          M2=MSFMIX(IJ,2)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            A12=(ASQ(IJ,1)-ASQ(IJ,2))*SINA*COSA 
     $                +ASQ(IJ,3)*(COSA2-SINA2)
            LAMB=SSXLAM(MH**2,M1**2,M2**2)         
            DWSF(IJ,3)=TEMP1*SQRT(LAMB)*A12**2
          ELSE IF(MH.LE.M12) THEN
               DWSF(IJ,3)=0.0
          ENDIF
C
          DWSF(IJ,4)=DWSF(IJ,3)
C
          IF(NUMH.EQ.1) THEN
            DO 121 JJ=1,4
              DWSFL(IJ,JJ)=DWSF(IJ,JJ)  
121         CONTINUE
          ELSE IF(NUMH.EQ.2) THEN
            DO 122 JJ=1,4
              DWSFH(IJ,JJ)=DWSF(IJ,JJ)  
122         CONTINUE
          ENDIF
C
150     CONTINUE
C
C          Now take care of sneutrinos. 
C    
        DO 155 IJ=10,12
          M1=MSFMIX(IJ,1)
          M2=MSFMIX(IJ,1)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSF(IJ,1)=SQRT(LAMB)*(ASQ(10,1))**2
     $                           /(16.0*PI*MH**3) 
          ELSE IF(MH.LE.M12) THEN
            DWSF(IJ,1) = 0.0
          ENDIF
          DWSF(IJ,2)=0.0
          DWSF(IJ,3)=0.0
          DWSF(IJ,4)=0.0
          IF(NUMH.EQ.1) THEN
            DO 151 JJ=1,4
              DWSFL(IJ,JJ)=DWSF(IJ,JJ)  
151         CONTINUE
          ELSE IF(NUMH.EQ.2) THEN
            DO 152 JJ=1,4
              DWSFH(IJ,JJ)=DWSF(IJ,JJ)  
152         CONTINUE
          ENDIF
C
155     CONTINUE 
        GO TO 1000
C
C
C          Pseudocalar neutral Higgses --> sfermions 
C           partial decay widths
C
233     TEMP1=GG/(2.0*AMW)
        BSQ(1)=TEMP1*AME*(EP-TANB*ATERM(1)) 
        BSQ(2)=TEMP1*AMMU*(EP-TANB*ATERM(2)) 
        BSQ(3)=TEMP1*AMTAU*(EP-TANB*ATERM(3)) 
        BSQ(4)=TEMP1*AMUP*(EP-COTB*ATERM(4)) 
        BSQ(5)=TEMP1*AMCH*(EP-COTB*ATERM(5)) 
        BSQ(6)=TEMP1*AMDN*(EP-TANB*ATERM(6)) 
        BSQ(7)=TEMP1*AMST*(EP-TANB*ATERM(7)) 
        BSQ(8)=TEMP1*MBQ*(EP-TANB*ATERM(8)) 
        BSQ(9)=TEMP1*MTQ*(EP-COTB*ATERM(9))
C
        DO 260 IJ=1,9
          IF(IJ.LT.4) THEN
            TEMP1=1.0/(16.0*PI*MH**3)
          ELSE 
            TEMP1=3.0/(16.0*PI*MH**3)
          ENDIF
          SINA=SIN(THETSF(IJ))
          COSA=COS(THETSF(IJ))
          SINA2=SINA**2
          COSA2=COSA**2
          M1=MSFMIX(IJ,1)
          M2=MSFMIX(IJ,1)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            B11=-2.0*COSA*SINA*BSQ(IJ)
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFP(IJ,1)=TEMP1*SQRT(LAMB)*B11**2
          ELSE IF(MH.LE.M12) THEN
            DWSFP(IJ,1)=0.0
          ENDIF
C
          M1=MSFMIX(IJ,2)
          M2=MSFMIX(IJ,2)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            B22=-B11
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFP(IJ,2)=TEMP1*SQRT(LAMB)*B22**2
          ELSE IF(MH.LE.M12) THEN
            DWSFP(IJ,2)=0.0
          ENDIF  
          M1=MSFMIX(IJ,1)
          M2=MSFMIX(IJ,2)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            B12=(COSA2-SINA2)*BSQ(IJ)
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFP(IJ,3)=TEMP1*SQRT(LAMB)*B12**2
          ELSE IF(MH.LE.M12) THEN
            DWSFP(IJ,3)=0.0
          ENDIF  
          DWSFP(IJ,4)=DWSFP(IJ,3)
260     CONTINUE
        DO 265 IJ=10,12
           DO 264 JJ=1,4
             DWSFP(IJ,JJ)=0.0
264        CONTINUE
265     CONTINUE
        GO TO 1000
C
C          Charged Higgses --> sfermions 
C           partial decay widths
C
333     TEMP1=-AMW*SIN2B
        CSQ(1,1)=GG*(TEMP1+(TANB*AMDN**2 + COTB*AMUP**2)/AMW)/SR2
        CSQ(2,1)=GG*(TEMP1+(TANB*AMST**2 + COTB*AMCH**2)/AMW)/SR2
        CSQ(3,1)=GG*(TEMP1+(TANB*MBQ**2 + COTB*MTQ**2)/AMW)/SR2
        CSQ(4,1)=GG*(TEMP1 + (TANB*AME**2)/AMW)/SR2
        CSQ(5,1)=GG*(TEMP1 + (TANB*AMMU**2)/AMW)/SR2
        CSQ(6,1)=GG*(TEMP1 + (TANB*AMTAU**2)/AMW)/SR2
C
        TEMP1=GG*(COTB+TANB)/(SR2*AMW)
        CSQ(1,2)=TEMP1*AMUP*AMDN         
        CSQ(2,2)=TEMP1*AMCH*AMST
        CSQ(3,2)=TEMP1*MTQ*MBQ
        CSQ(4,2)=0.0
        CSQ(5,2)=0.0
        CSQ(6,2)=0.0
C
        TEMP1=GG/(SR2*AMW)
        CSQ(1,3)=TEMP1*AMUP*(EP-COTB*ATERM(4))
        CSQ(2,3)=TEMP1*AMCH*(EP-COTB*ATERM(5))
        CSQ(3,3)=TEMP1*MTQ*(EP-COTB*ATERM(9))
        CSQ(4,3)=0.0
        CSQ(5,3)=0.0
        CSQ(6,3)=0.0
C
        CSQ(1,4)=TEMP1* AMDN*(EP-TANB*ATERM(6))
        CSQ(2,4)=TEMP1* AMST*(EP-TANB*ATERM(7))
        CSQ(3,4)=TEMP1* MBQ*(EP-TANB*ATERM(8))
        CSQ(4,4)=TEMP1* AME*(EP-TANB*ATERM(1))
        CSQ(5,4)=TEMP1* AMMU*(EP-TANB*ATERM(2))
        CSQ(6,4)=TEMP1* AMTAU*(EP-TANB*ATERM(3))
C
        DO 350 IC=1,3
          TEMP1=3.0/(16.0*PI*MH**3)
          IF(IC.EQ.1) THEN
            IJU=4
            IJD=6
          ELSE IF(IC.EQ.2) THEN
            IJU=5
            IJD=7
          ELSE IF(IC.EQ.3) THEN
            IJU=9
            IJD=8
          ENDIF
          SINAU=SIN(THETSF(IJU))
          COSAU=COS(THETSF(IJU))
          SINAD=SIN(THETSF(IJD))
          COSAD=COS(THETSF(IJD))
C
          M1=MSFMIX(IJU,1)
          M2=MSFMIX(IJD,1)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            C11=COSAU*COSAD*CSQ(IC,1) 
     $            + SINAU*SINAD*CSQ(IC,2)
     $                - SINAU*COSAD*CSQ(IC,3) 
     $                    - COSAU*SINAD*CSQ(IC,4)     
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFC(IC,1)=TEMP1*SQRT(LAMB)*C11**2
          ELSE IF(MH.LE.M12) THEN
            DWSFC(IC,1) = 0.0
          ENDIF  
C
          M1=MSFMIX(IJU,1)
          M2=MSFMIX(IJD,2)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            C12=COSAU*SINAD*CSQ(IC,1) 
     $            - SINAU*COSAD*CSQ(IC,2)
     $                - SINAU*SINAD*CSQ(IC,3) 
     $                    + COSAU*COSAD*CSQ(IC,4)
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFC(IC,2)=TEMP1*SQRT(LAMB)*C12**2
          ELSE IF(MH.LE.M12) THEN
            DWSFC(IC,2)=0.0
          ENDIF  
C
          M1=MSFMIX(IJU,2)
          M2=MSFMIX(IJD,1)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            C21=SINAU*COSAD*CSQ(IC,1) 
     $            - COSAU*SINAD*CSQ(IC,2)
     $                + COSAU*COSAD*CSQ(IC,3) 
     $                    - SINAU*SINAD*CSQ(IC,4) 
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFC(IC,3)=TEMP1*SQRT(LAMB)*C21**2
          ELSE IF(MH.LE.M12) THEN
            DWSFC(IC,3)=0.0
          ENDIF  
C
          M1=MSFMIX(IJU,2)
          M2=MSFMIX(IJD,2)
          M12=M1+M2
          IF(MH.GT.M12) THEN
            C22=SINAU*SINAD*CSQ(IC,1) 
     $            + COSAU*COSAD*CSQ(IC,2)
     $                + COSAU*SINAD*CSQ(IC,3) 
     $                    - SINAU*COSAD*CSQ(IC,4) 
            LAMB=SSXLAM(MH**2,M1**2,M2**2)
            DWSFC(IC,4)=TEMP1*SQRT(LAMB)*C22**2
          ELSE IF(MH.LE.M12) THEN
            DWSFC(IC,4)=0.0
          ENDIF
C
350     CONTINUE
C
C
C         Now calculate the sleptonic
C          partial decay widths of the 
C          charged Higgs.
C
        DO 355 IC = 4,6
          TEMP1=1.0/(16.0*PI*MH**3)
            IF(IC.EQ.4) THEN
              IJU=10
              IJD=1
            ELSE IF(IC.EQ.5) THEN
              IJU=11
              IJD=2
            ELSE IF(IC.EQ.6) THEN
              IJU=12
              IJD=3
            ENDIF
            SINAD=SIN(THETSF(IJD))
            COSAD=COS(THETSF(IJD))
C
            M1=MSFMIX(IJU,1)
            M2=MSFMIX(IJD,1)
            M12=M1+M2
            IF(MH.GT.M12) THEN
              C11=COSAD*CSQ(IC,1)-SINAD*CSQ(IC,4)
              LAMB=SSXLAM(MH**2,M1**2,M2**2)
              DWSFC(IC,1)=TEMP1*SQRT(LAMB)*C11**2
            ELSE IF(MH.LE.M12) THEN
              DWSFC(IC,1)=0.0
            ENDIF
C
            M1=MSFMIX(IJU,1)
            M2=MSFMIX(IJD,2)
            M12=M1+M2
            IF(MH.GT.M12) THEN
              C12=SINAD*CSQ(IC,1)+COSAD*CSQ(IC,4)  
              LAMB=SSXLAM(MH**2,M1**2,M2**2)
              DWSFC(IC,2)=TEMP1*SQRT(LAMB)*C12**2
            ELSE IF(MH.LE.M12) THEN
              DWSFC(IC,2)=0.0
            ENDIF  
            DWSFC(IC,3)=0.0
            DWSFC(IC,4)=0.0
355     CONTINUE
1000  CONTINUE
C          H_l decays
      CALL SSSAVE(ISHL,DWSFL(1,1),ISEL,-ISEL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(1,2),ISER,-ISER,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(2,1),ISMUL,-ISMUL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(2,2),ISMUR,-ISMUR,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(3,1),ISTAU1,-ISTAU1,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(3,2),ISTAU2,-ISTAU2,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(3,3),ISTAU1,-ISTAU2,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(3,4),ISTAU2,-ISTAU1,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(4,1),ISUPL,-ISUPL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(4,2),ISUPR,-ISUPR,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(5,1),ISCHL,-ISCHL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(5,2),ISCHR,-ISCHR,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(6,1),ISDNL,-ISDNL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(6,2),ISDNR,-ISDNR,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(7,1),ISSTL,-ISSTL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(7,2),ISSTR,-ISSTR,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(8,1),ISBT1,-ISBT1,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(8,2),ISBT2,-ISBT2,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(8,3),ISBT1,-ISBT2,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(8,4),ISBT2,-ISBT1,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(9,1),ISTP1,-ISTP1,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(9,2),ISTP2,-ISTP2,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(9,3),ISTP1,-ISTP2,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(9,4),ISTP2,-ISTP1,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(10,1),ISNEL,-ISNEL,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(11,1),ISNML,-ISNML,0,0,0)
      CALL SSSAVE(ISHL,DWSFL(12,1),ISNTL,-ISNTL,0,0,0)
C         H_h decays
      CALL SSSAVE(ISHH,DWSFH(1,1),ISEL,-ISEL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(1,2),ISER,-ISER,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(2,1),ISMUL,-ISMUL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(2,2),ISMUR,-ISMUR,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(3,1),ISTAU1,-ISTAU1,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(3,2),ISTAU2,-ISTAU2,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(3,3),ISTAU1,-ISTAU2,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(3,4),ISTAU2,-ISTAU1,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(4,1),ISUPL,-ISUPL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(4,2),ISUPR,-ISUPR,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(5,1),ISCHL,-ISCHL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(5,2),ISCHR,-ISCHR,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(6,1),ISDNL,-ISDNL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(6,2),ISDNR,-ISDNR,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(7,1),ISSTL,-ISSTL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(7,2),ISSTR,-ISSTR,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(8,1),ISBT1,-ISBT1,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(8,2),ISBT2,-ISBT2,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(8,3),ISBT1,-ISBT2,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(8,4),ISBT2,-ISBT1,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(9,1),ISTP1,-ISTP1,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(9,2),ISTP2,-ISTP2,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(9,3),ISTP1,-ISTP2,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(9,4),ISTP2,-ISTP1,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(10,1),ISNEL,-ISNEL,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(11,1),ISNML,-ISNML,0,0,0)
      CALL SSSAVE(ISHH,DWSFH(12,1),ISNTL,-ISNTL,0,0,0)
C          Decay of H_p
      CALL SSSAVE(ISHA,DWSFP(1,3),ISEL,-ISER,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(1,4),ISER,-ISEL,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(2,3),ISMUL,-ISMUR,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(2,4),ISMUR,-ISMUL,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(3,1),ISTAU1,-ISTAU1,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(3,2),ISTAU2,-ISTAU2,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(3,3),ISTAU1,-ISTAU2,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(3,4),ISTAU2,-ISTAU1,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(4,3),ISUPL,-ISUPR,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(4,4),ISUPR,-ISUPL,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(5,3),ISCHL,-ISCHR,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(5,4),ISCHR,-ISCHL,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(6,3),ISDNL,-ISDNR,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(6,4),ISDNR,-ISDNL,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(7,3),ISSTL,-ISSTR,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(7,4),ISSTR,-ISSTL,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(8,1),ISBT1,-ISBT1,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(8,2),ISBT2,-ISBT2,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(8,3),ISBT1,-ISBT2,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(8,4),ISBT2,-ISBT1,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(9,1),ISTP1,-ISTP1,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(9,2),ISTP2,-ISTP2,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(9,3),ISTP1,-ISTP2,0,0,0)
      CALL SSSAVE(ISHA,DWSFP(9,4),ISTP2,-ISTP1,0,0,0)
C          Decay of H+
      CALL SSSAVE(ISHC,DWSFC(1,1),ISUPL,-ISDNL,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(1,2),ISUPR,-ISDNR,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(2,1),ISCHL,-ISSTL,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(2,2),ISCHR,-ISSTR,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(3,1),ISTP1,-ISBT1,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(3,2),ISTP1,-ISBT2,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(3,3),ISTP2,-ISBT1,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(3,4),ISTP2,-ISBT2,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(4,1),-ISEL,ISNEL,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(5,1),-ISMUL,ISNML,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(6,1),-ISTAU1,ISNTL,0,0,0)
      CALL SSSAVE(ISHC,DWSFC(6,2),-ISTAU2,ISNTL,0,0,0)
      RETURN
      END

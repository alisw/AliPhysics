#include "isajet/pilot.h"
      FUNCTION LABEL(ID)
C
C          Return the CHARACTER*8 label for the particle ID.
C          Quark-based IDENT code.
C          MSSM names for squarks, sleptons, Higgs bosons.
C
C          Ver. 7.49: Offset of INDEX must match that in FLAVOR.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/qlmass.inc"
C
#if defined(CERNLIB_LEVEL2)
C          DUMMY COMMON BLOCK TO ALLOW LEVEL2 STORAGE.
      COMMON/XLABEL/LLEP,LMES0,LMES1,LBAR0,LABAR0,LBAR1,LABAR1,LQQ,LAQQ
      LEVEL2,/XLABEL/
#endif
      INTEGER ID
      CHARACTER*8 LABEL
      CHARACTER*8 LLEP,LMES0,LMES1,LBAR0,LABAR0,LBAR1,LABAR1
      CHARACTER*8 LQQ,LAQQ
      DIMENSION LLEP(149)
      DIMENSION LMES0(64),LMES1(64)
      DIMENSION LBAR0(109),LABAR0(109),LBAR1(109),LABAR1(109)
      DIMENSION LQQ(21),LAQQ(21)
      INTEGER IFL1,IFL2,IFL3,JSPIN,INDEX,I,J,IDABS
C
C          Diquark labels
C
      DATA LQQ/
     1'UU0. ','UD0. ','DD0. ','US0. ','DS0. ','SS0. ','UC0. ','DC0. ',
     2'SC0. ','CC0. ','UB0. ','DB0. ','SB0. ','CB0. ','BB0. ','UT0. ',
     3'DT0. ','ST0. ','CT0. ','BT0. ','TT0. '/
      DATA LAQQ/
     1'AUU0.','AUD0.','ADD0.','AUS0.','ADS0.','ASS0.','AUC0.','ADC0.',
     2'ASC0.','ACC0.','AUB0.','ADB0.','ASB0.','ACB0.','ABB0.','AUT0.',
     3'ADT0.','AST0.','ACT0.','ABT0.','ATT0.'/
C
C          Quark and lepton labels
C
      DATA LLEP/
     $'     ','UP   ','UB   ','DN   ','DB   ','ST   ','SB   ','CH   ',
     $'CB   ','BT   ','BB   ','TP   ','TB   ','Y    ','YB   ','X    ',
     $'XB   ','GL   ','ERR  ','GM   ','ERR  ','NUE  ','ANUE ','E-   ',
     $'E+   ','NUM  ','ANUM ','MU-  ','MU+  ','NUT  ','ANUT ','TAU- ',
     $'TAU+ ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','KS   ',
     $'ERR  ','ERR  ','KL   ',
     $'UPL  ','UBL  ','DNL  ','DBL  ','STL  ','SBL  ','CHL  ','CBL  ',
     $'BT1  ','BB1  ','TP1  ','TB1  ','ERR  ','ERR  ','ERR  ','ERR  ',
     $'GLSS ','ERR  ','Z1SS ','ERR  ','NUEL ','ANUEL','EL-  ','EL+  ',
     $'NUML ','ANUML','MUL- ','MUL+ ','NUTL ','ANUTL','TAU1-','TAU1+',
     $'ERR  ','ERR  ','ERR  ','ERR  ','W1SS+','W1SS-','Z2SS ','ERR  ',
     $'UPR  ','UBR  ','DNR  ','DBR  ','STR  ','SBR  ','CHR  ','CBR  ',
     $'BT2  ','BB2  ','TP2  ','TB2  ','ERR  ','ERR  ','ERR  ','ERR  ',
     $'W2SS+','W2SS-','Z3SS ','ERR  ','NUER ','ANUER','ER-  ','ER+  ',
     $'NUMR ','ANUMR','MUR- ','MUR+ ','NUTR ','ANUTR','TAU2-','TAU2+',
     $'ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','Z4SS ','ERR  ',
     $'W+   ','W-   ','HIGGS','ERR  ','HL0  ','ERR  ','HH0  ','ERR  ',
     $'HA0  ','ERR  ','H40  ','AH40 ','H+   ','H-   ','H2+  ','H2-  ',
     $'H1++ ','H1-- ','H2++ ','H2-- ','Z0   ','ERR  ','GVSS ','ERR  ',
     $'GRAV ','ERR  '/
C
C          0- meson labels
C
      DATA LMES0/
     1'PI0  ','PI+  ','ETA  ','PI-  ','K+   ','K0   ','ETAP ','AK0  ',
     2'K-   ','AD0  ','D-   ','DS-  ','ETAC ','DS+  ','D+   ','D0   ',
     2'B+   ','B0   ','BS   ','BC   ','ETAB ','ABC  ','ABS  ','AB0  ',
     3'B-   ','UT.  ','DT.  ','ST.  ','CT.  ','BT.  ','TT.  ','TB.  ',
     4'TC.  ','TS.  ','TD.  ','TU.  ','UY.  ','DY.  ','SY.  ','CY.  ',
     5'BY.  ','TY.  ','YY.  ','YT.  ','YB.  ','YC.  ','YS.  ','YD.  ',
     6'YU.  ','UX.  ','DX.  ','SX.  ','CX.  ','BX.  ','TX.  ','YX.  ',
     7'XX.  ','XY.  ','XT.  ','XB.  ','XC.  ','XS.  ','XD.  ','XU.  '/
C
C          1- meson labels
C
      DATA LMES1/
     1'RHO0 ','RHO+ ','OMEG ','RHO- ','K*+  ','K*0  ','PHI  ','AK*0 ',
     2'K*-  ','AD*0 ','D*-  ','DS*- ','JPSI ','DS*+ ','D*+  ','D*0  ',
     3'B*+  ','B*0  ','BS*  ','BC*  ','UPSL ','ABC* ','ABS* ','AB0* ',
     4'B*-  ','UT*  ','DT*  ','ST*  ','CT*  ','BT*  ','TT*  ','TB*  ',
     5'TC*  ','TS*  ','TD*  ','TU*  ','UY*  ','DY*  ','SY*  ','CY*  ',
     6'BY*  ','TY*  ','YY*  ','YT*  ','YB*  ','YC*  ','YS*  ','YD*  ',
     7'YU*  ','UX*  ','DX*  ','SX*  ','CX*  ','BX*  ','TX*  ','YX*  ',
     8'XX*  ','XY*  ','XT*  ','XB*  ','XC*  ','XS*  ','XD*  ','XU*  '/
C
C          1/2+ baryon labels
C
      DATA LBAR0/
     1'ERR  ','P    ','N    ','ERR  ','ERR  ','S+   ','S0   ','S-   ',
     2'L    ','XI0  ','XI-  ','ERR  ','ERR  ','ERR  ','SC++ ','SC+  ',
     3'SC0  ','LC+  ','USC. ','DSC. ','SSC. ','SDC. ','SUC. ','UCC. ',
     4'DCC. ','SCC. ','ERR  ','ERR  ','ERR  ','ERR  ','UUB. ','UDB. ',
     5'DDB. ','DUB. ','USB. ','DSB. ','SSB. ','SDB. ','SUB. ','UCB. ',
     6'DCB. ','SCB. ','CCB. ','CSB. ','CDB. ','CUB. ','UBB. ','DBB. ',
     7'SBB. ','CBB. ','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','UUT. ',
     8'UDT. ','DDT. ','DUT. ','UST. ','DST. ','SST. ','SDT. ','SUT. ',
     9'UCT. ','DCT. ','SCT. ','CCT. ','CST. ','CDT. ','CUT. ','UBT. ',
     1'DBT. ','SBT. ','CBT. ','BBT. ','BCT. ','BST. ','BDT. ','BUT. ',
     2'UTT. ','DTT. ','STT. ','CTT. ','BTT. ','ERR  ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','UUY. ','UDY. ','DDY. ','DUY. ','USY. ',
     4'DSY. ','SSY. ','SDY. ','SUY. ','UUX. ','UDX. ','DDX. ','DUX. ',
     5'USX. ','DSX. ','SSX. ','SDX. ','SUX. '/
      DATA LABAR0/
     1'ERR  ','AP   ','AN   ','ERR  ','ERR  ','AS-  ','AS0  ','AS+  ',
     2'AL   ','AXI0 ','AXI+ ','ERR  ','ERR  ','ERR  ','ASC--','ASC- ',
     3'ASC0 ','ALC- ','AUSC.','ADSC.','ASSC.','ASDC.','ASUC.','AUCC.',
     4'ADCC.','ASCC.','ERR  ','ERR  ','ERR  ','ERR  ','AUUB.','AUDB.',
     5'ADDB.','ADUB.','AUSB.','ADSB.','ASSB.','ASDB.','ASUB.','AUCB.',
     6'ADCB.','ASCB.','ACCB.','ACSB.','ACDB.','ACUB.','AUBB.','ADBB.',
     7'ASBB.','ACBB.','ERR  ','ERR  ','ERR  ','ERR  ','ERR  ','AUUT.',
     8'AUDT.','ADDT.','ADUT.','AUST.','ADST.','ASST.','ASDT.','ASUT.',
     9'AUCT.','ADCT.','ASCT.','ACCT.','ACST.','ACDT.','ACUT.','AUBT.',
     1'ADBT.','ASBT.','ACBT.','ABBT.','ABCT.','ABST.','ABDT.','ABUT.',
     2'AUTT.','ADTT.','ASTT.','ACTT.','ABTT.','ERR  ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','AUUY.','AUDY.','ADDY.','ADUY.','AUSY.',
     4'ADSY.','ASSY.','ASDY.','ASUY.','AUUX.','AUDX.','ADDX.','ADUX.',
     5'AUSX.','ADSX.','ASSX.','ASDX.','ASUX.'/
C
C          3/2+ baryon labels
C
      DATA LBAR1/
     1'DL++ ','DL+  ','DL0  ','DL-  ','ERR  ','S*+  ','S*0  ','S*-  ',
     2'ERR  ','XI*0 ','XI*- ','OM-  ','ERR  ','ERR  ','UUC* ','UDC* ',
     3'DDC* ','ERR  ','USC* ','DSC* ','SSC* ','ERR  ','ERR  ','UCC* ',
     4'DCC* ','SCC* ','CCC* ','ERR  ','ERR  ','ERR  ','UUB* ','UDB* ',
     5'DDB* ','ERR  ','USB* ','DSB* ','SSB* ','ERR  ','ERR  ','UCB* ',
     6'DCB* ','SCB* ','CCB* ','ERR  ','ERR  ','ERR  ','UBB* ','DBB* ',
     7'SBB* ','CBB* ','BBB* ','ERR  ','ERR  ','ERR  ','ERR  ','UUT* ',
     8'UDT* ','DDT* ','ERR  ','UST* ','DST* ','SST* ','ERR  ','ERR  ',
     9'UCT* ','DCT* ','SCT* ','CCT* ','ERR  ','ERR  ','ERR  ','UBT* ',
     1'DBT* ','SBT* ','CBT* ','BBT* ','ERR  ','ERR  ','ERR  ','ERR  ',
     2'UTT* ','DTT* ','STT* ','CTT* ','BTT* ','TTT* ','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','UUY* ','UDY* ','DDY* ','ERR  ','USY* ',
     4'DSY* ','SSY* ','ERR  ','ERR  ','UUX* ','UDX* ','DDX* ','ERR  ',
     5'USX* ','DSX* ','SSX* ','ERR  ','ERR  '/
      DATA LABAR1/
     1'ADL--','ADL- ','ADL0 ','ADL+ ','ERR  ','AS*- ','AS*0 ','AS*+ ',
     2'ERR  ','AXI*0','AXI*+','AOM+ ','ERR  ','ERR  ','AUUC*','AUDC*',
     3'ADDC*','ERR  ','AUSC*','ADSC*','ASSC*','ERR  ','ERR  ','AUCC*',
     4'ADCC*','ASCC*','ACCC*','ERR  ','ERR  ','ERR  ','AUUB*','AUDB*',
     5'ADDB*','ERR  ','AUSB*','ADSB*','ASSB*','ERR  ','ERR  ','AUCB*',
     6'ADCB*','ASCB*','ACCB*','ERR  ','ERR  ','ERR  ','AUBB*','ADBB*',
     7'ASBB*','ACBB*','ABBB*','ERR  ','ERR  ','ERR  ','ERR  ','AUUT*',
     8'AUDT*','ADDT*','ERR  ','AUST*','ADST*','ASST*','ERR  ','ERR  ',
     9'AUCT*','ADCT*','ASCT*','ACCT*','ERR  ','ERR  ','ERR  ','AUBT*',
     1'ADBT*','ASBT*','ACBT*','ABBT*','ERR  ','ERR  ','ERR  ','ERR  ',
     2'AUTT*','ADTT*','ASTT*','ACTT*','ABTT*','ATTT*','ERR  ','ERR  ',
     3'ERR  ','ERR  ','ERR  ','AUUY*','AUDY*','ADDY*','ERR  ','AUSY*',
     4'ADSY*','ASSY*','ERR  ','ERR  ','AUUX*','AUDX*','ADDX*','ERR  ',
     5'AUSX*','ADSX*','ASSX*','ERR  ','ERR  '/
C
C          Entry
C
      LABEL='ERR'
      IDABS=IABS(ID)
      IF(IDABS.EQ.0) THEN
        LABEL='     '
        RETURN
      ENDIF
      CALL FLAVOR(ID,IFL1,IFL2,IFL3,JSPIN,INDEX)
      IF(INDEX.LE.0) RETURN
      IF(IDABS.GT.10000.OR.JSPIN.GT.1) GO TO 500
      IF(IDABS.LT.100) GO TO 200
      IF(IDABS.LT.1000) GO TO 100
      IF(ID.NE.0.AND.MOD(ID,100).EQ.0) GO TO 300
C
C          Baryons
C
      INDEX=INDEX-109*JSPIN-36*NMES-NQLEP
      INDEX=INDEX-13
      IF(JSPIN.EQ.0.AND.ID.GT.0) LABEL=LBAR0(INDEX)
      IF(JSPIN.EQ.0.AND.ID.LT.0) LABEL=LABAR0(INDEX)
      IF(JSPIN.EQ.1.AND.ID.GT.0) LABEL=LBAR1(INDEX)
      IF(JSPIN.EQ.1.AND.ID.LT.0) LABEL=LABAR1(INDEX)
      GO TO 999
C
C          Mesons
C
100   CONTINUE
      I=MAX0(IFL2,IFL3)
      J=-MIN0(IFL2,IFL3)
      INDEX=MAX0(I-1,J-1)**2+I+MAX0(I-J,0)
      IF(JSPIN.EQ.0) LABEL=LMES0(INDEX)
      IF(JSPIN.EQ.1) LABEL=LMES1(INDEX)
      GO TO 999
C
C          Quarks, leptons, etc.
C
200   CONTINUE
      INDEX=2*INDEX
      IF(ID.LE.0) INDEX=INDEX+1
      LABEL=LLEP(INDEX)
      GO TO 999
300   I=IABS(IFL1)
      J=IABS(IFL2)
      INDEX=I+J*(J-1)/2
      IF(ID.GT.0) LABEL=LQQ(INDEX)
      IF(ID.LT.0) LABEL=LAQQ(INDEX)
      RETURN
C
C          Special hadrons - used only in B decays
C
500   CONTINUE
      IF(ID.EQ.10121) THEN
        LABEL='A1+'
      ELSEIF(ID.EQ.-10121) THEN
        LABEL='A1-'
      ELSEIF(ID.EQ.10111) THEN
        LABEL='A10'
      ELSEIF(ID.EQ.10131) THEN
        LABEL='K1+'
      ELSEIF(ID.EQ.-10131) THEN
        LABEL='K1-'
      ELSEIF(ID.EQ.10231) THEN
        LABEL='K10'
      ELSEIF(ID.EQ.-10231) THEN
        LABEL='AK10'
      ELSEIF(ID.EQ.30131) THEN
        LABEL='K1*+'
      ELSEIF(ID.EQ.-30131) THEN
        LABEL='K1*-'
      ELSEIF(ID.EQ.30231) THEN
        LABEL='K1*0'
      ELSEIF(ID.EQ.-30231) THEN
        LABEL='AK1*0'
      ELSEIF(ID.EQ.132) THEN
        LABEL='K2*+'
      ELSEIF(ID.EQ.-132) THEN
        LABEL='K2*-'
      ELSEIF(ID.EQ.232) THEN
        LABEL='K2*0'
      ELSEIF(ID.EQ.-232) THEN
        LABEL='AK2*0'
      ELSEIF(ID.EQ.10110) THEN
        LABEL='F0'
      ELSEIF(ID.EQ.112) THEN
        LABEL='F2'
      ELSEIF(ID.EQ.10441) THEN
        LABEL='PSI2'
      ELSEIF(ID.EQ.20440) THEN
        LABEL='CHI0'
      ELSEIF(ID.EQ.20441) THEN
        LABEL='CHI1'
      ELSEIF(ID.EQ.20442) THEN
        LABEL='CHI2'
      ELSE
        LABEL='ERR'
      ENDIF
999   RETURN
      END

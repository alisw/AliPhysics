*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:14  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      SUBROUTINE PCSDAT(LUNIN,LUNOUT,PRFLAG)
C
C *** PREPARATION OF DATA STMTS. FOR ELASTIC AND INELASTIC MEASURED  ***
C *** CROSS SECTIONS OF PION,KAON AND PROTON/ANTIPROTON/NEUTRON      ***
C *** ON PROTONS. CALCULATE FROM THIS CROSS SECTIONS FOR STRANGE     ***
C *** BARYONS ON PROTONS USING ADDITIVE QUARK-QUARK SCATTERING MODEL ***
C *** NVE 22-FEB-1988 CERN GENEVA ***
C
C THE PROGRAM PRODUCES AN OUTPUT FILE WHICH CONTAINS THE DATA
C STATEMENTS FOR ALL THE NEEDED CROSS-SECTIONS IN SEQUENCE "/CSDAT"
C THESE DATA STATEMENTS MAY BE USED FOR ROUTINE "GHESIG"
C
C LUNIN  = UNIT NUMBER FOR CROSS-SECTION DATA FILE
C LUNOUT = UNIT NUMBER FOR DATA STATEMENTS TO BE WRITTEN
C PRFLAG = PRINTOUT FLAG (TRUE/FALSE)
C
C ORIGIN : H.FESEFELDT 06-OCT-87 (ROUTINE CSDATA)
C
C QUARK SCATTERING AMPLITUDES
C <PP |PP > = <NN |NN > = P
C <PN |PN > = <NP |NP > = P
C <PN |NP > = <NP |PN > = 0
C <PL |PL > = <NL |NL > = P - S
C <PLB|PLB> = <NLB|NLB> = P - S
C <PL |LP > = <NL |LN > = 0
C <PBN|PBN> = <NBP|NBP> = P
C <PBP|PBP> = <NBN|NBN> = P + A
C <PBP|NBN> = <NBN|PBP> = A
C <LL |LL > =           = (P - S)**2/P
C <LBL|LBL> =           = (P + A)*(P + S)**2/P**2
C <PBP|LBL> = <NBN|LBL> = A*(P - S)/P
C
C CROSS SECTIONS
C PI-  P  = 6P + 2A
C PI+  P  = 6P +  A
C K-   P  = 6P + 2A - 3S
C K0S  P  = 6P + A/2- 3S
C K0L  P  = 6P + A/2- 3S
C K0   P  = 6P      - 3S  ==> S(K0)=S(K+)
C K0B  P  = 6P + A  - 3S  ==> S(K0B)=S(K+)/4+S(K0L)/3+5S(K-)/12
C K+   P  = 6P      - 3S
C P    P  = 9P
C PB   P  = 9P + 5A
C N    P  = 9P
C
C FOR THE FOLLOWING REACTIONS WE TOOK THE AMPLITUDES
C
C P = S(P P)/9
C A = (S(PB P)-S(P P))/5
C S = 2S(P P)/9-S(K+ P)/3
C
C =====>
C
C NB   P  = 9P + 4A
C L    P  = 9P      - 3S
C LB   P  = 9P + 2A - 3S
C S-   P  = 9P      - 3S
C S+   P  = 9P      - 3S
C S-B  P  = 9P + 2A - 3S
C S+B  P  = 9P + 4A - 3S
C X0   P  = 9P      - 6S
C X-   P  = 9P      - 6S
C X0B  P  = 9P + 2A - 6S
C X-B  P  = 9P +  A - 6S
C
#include "geant321/s_csdim.inc"
C
      LOGICAL PRFLAG
C
      DIMENSION A(6,41)
      CHARACTER*72 IA
      DIMENSION AX(9),LX(3),MX(3),BX(3),CX(9),EX(9)
      DATA AX/423.,78.,-54.,78.,34.25,-7.5,-54.,-7.5,27./
C
C --- INITIALIZE ALL ARRAYS ---
      DO 8000 J=1,41
      PLAB(J)=0.0
      DO 8001 I=1,35
      CSEL(I,J)=0.0
      CSIN(I,J)=0.0
      IF (I .GT. 3) GO TO 8001
      CSPIEL(I,J)=0.0
      CSPIIN(I,J)=0.0
      CSPNEL(I,J)=0.0
      CSPNIN(I,J)=0.0
 8001 CONTINUE
 8000 CONTINUE
C
      DO 8002 J=1,17
      ELAB(J)=0.0
      DO 8003 I=1,15
      CNLWEL(I,J)=0.0
      CNLWIN(I,J)=0.0
 8003 CONTINUE
      IF (J .GT. 15) GO TO 8002
      CNLWAT(J)=0.0
 8002 CONTINUE
C
      DO 8004 J=1,21
      EKFISS(J)=0.0
      DO 8005 I=1,4
      CSFISS(I,J)=0.0
 8005 CONTINUE
 8004 CONTINUE
C
      DO 8006 I=1,100
      CSCAP(I)=0.0
 8006 CONTINUE
C
      IP=0
      IAT=-3
      CNLWAT( 1)=   1.
      CNLWAT( 2)=  16.
      CNLWAT( 3)=  27.
      CNLWAT( 4)=  56.
      CNLWAT( 5)=  59.
      CNLWAT( 6)=  64.
      CNLWAT( 7)=  91.
      CNLWAT( 8)= 112.
      CNLWAT( 9)= 119.
      CNLWAT(10)= 127.
      CNLWAT(11)= 137.
      CNLWAT(12)= 181.
      CNLWAT(13)= 207.
      CNLWAT(14)= 209.
      CNLWAT(15)= 238.
      CALL UCOPY(AX,CX,9)
      CALL MINV(AX,3,DET,LX,MX)
      EX(1)=AX(1)*CX(1)+AX(4)*CX(2)+AX(7)*CX(3)
      EX(2)=AX(1)*CX(4)+AX(4)*CX(5)+AX(7)*CX(6)
      EX(3)=AX(1)*CX(7)+AX(4)*CX(8)+AX(7)*CX(9)
      EX(4)=AX(2)*CX(1)+AX(5)*CX(2)+AX(8)*CX(3)
      EX(5)=AX(2)*CX(4)+AX(5)*CX(5)+AX(8)*CX(6)
      EX(6)=AX(2)*CX(7)+AX(5)*CX(8)+AX(8)*CX(9)
      EX(7)=AX(3)*CX(1)+AX(6)*CX(2)+AX(9)*CX(3)
      EX(8)=AX(3)*CX(4)+AX(6)*CX(5)+AX(9)*CX(6)
      EX(9)=AX(3)*CX(7)+AX(6)*CX(8)+AX(9)*CX(9)
      IF(PRFLAG) PRINT 7001,AX,DET,EX
 7001 FORMAT(1H ,'CSDATA MATRIX INVERSION ',9F10.6,2X,F10.1/
     *       1H ,'       UNIT MATRIX      ',9F10.6)
    1 READ(LUNIN,1001) IA
      IF (IA(1:3) .EQ. 'END') GO TO 100
      IF(PRFLAG) PRINT 1002,IA
      READ(LUNIN,1001) IA
      IF(PRFLAG) PRINT 1002,IA
      READ(LUNIN,1001) IA
      IF(PRFLAG) PRINT 1002,IA
 1001 FORMAT(A)
 1002 FORMAT(1H ,3X,A)
      IP=IP+1
      IF(IP.GT.11) GOTO 95
      IF(IP.GT.10) GOTO 90
      IF(IP.GT. 5) GOTO 3
      DO 2 I=1,41
      READ(LUNIN,1003) PLAB(I),(A(J,I),J=1,6)
      IF(PRFLAG) PRINT 1004,PLAB(I),(A(J,I),J=1,6)
 1003 FORMAT(7F10.2)
 1004 FORMAT(1H ,7F10.2)
    2 CONTINUE
      GOTO 5
    3 READ(LUNIN,1001) IA
      IF(PRFLAG) PRINT 1002,IA
      DO 4 I=1,17
      READ(LUNIN,1003) ELAB(I),(A(J,I),J=1,6)
      IF(PRFLAG) PRINT 1010,ELAB(I),(A(J,I),J=1,6)
 1010 FORMAT(1H ,F10.4,6F10.2)
    4 CONTINUE
    5 IF(IP.EQ. 1) GOTO 20
      IF(IP.EQ. 2) GOTO 30
      IF(IP.EQ. 3) GOTO 40
      IF(IP.EQ. 4) GOTO 60
      IF(IP.EQ. 5) GOTO 70
      IF(IP.GE. 6) GOTO 80
      GOTO 1
   20 DO 21 I=1,41
      CSEL(9,I)=A(1,I)
      CSIN(9,I)=A(2,I)
      CSEL(7,I)=A(3,I)
   21 CSIN(7,I)=A(4,I)
      GOTO 1
   30 DO 31 I=1,41
      CSEL(13,I)=A(1,I)
      CSIN(13,I)=A(2,I)
      CSEL(10,I)=A(3,I)
      CSIN(10,I)=A(4,I)
      CSEL(12,I)=A(5,I)
      CSIN(12,I)=A(6,I)
   31 CONTINUE
      GOTO 1
   40 IF(PRFLAG) PRINT 1008
 1008 FORMAT(1H0,'QUARK PARTON SCATTERING AMPLITUDES'/
     *1H ,'    P(GEV/C)     P         A         S          P          A
     *         S'/
     *1H ,'                 LEAST SQUARE FIT               CALCULATED')
      DO 41 I=1,41
      CSEL(14,I)=A(1,I)
      CSIN(14,I)=A(2,I)
      CSEL(15,I)=A(3,I)
      CSIN(15,I)=A(4,I)
      CSEL(16,I)=A(5,I)
      CSIN(16,I)=A(6,I)
      BX(1)=  6.*(CSEL( 7,I)+CSIN( 7,I))+ 6.*(CSEL( 9,I)+CSIN( 9,I))
     *      + 6.*(CSEL(10,I)+CSIN(10,I))+ 6.*(CSEL(12,I)+CSIN(12,I))
     *      + 6.*(CSEL(13,I)+CSIN(13,I))+ 9.*(CSEL(14,I)+CSIN(14,I))
     *      + 9.*(CSEL(15,I)+CSIN(15,I))+ 9.*(CSEL(16,I)+CSIN(16,I))
      BX(2)=     (CSEL( 7,I)+CSIN( 7,I))+    (CSEL( 9,I)+CSIN( 9,I))
     *      + .5*(CSEL(12,I)+CSIN(12,I))+ 5.*(CSEL(15,I)+CSIN(15,I))
      BX(3)=- 3.*(CSEL(10,I)+CSIN(10,I))- 3.*(CSEL(12,I)+CSIN(12,I))
     *      - 3.*(CSEL(13,I)+CSIN(13,I))
      PAM =AX(1)*BX(1)+AX(4)*BX(2)+AX(7)*BX(3)
      AAM =AX(2)*BX(1)+AX(5)*BX(2)+AX(8)*BX(3)
      SAM =AX(3)*BX(1)+AX(6)*BX(2)+AX(9)*BX(3)
C     RATKP=CSEL(10,I)/(CSEL(10,I)+CSIN(10,I))
C     RATKM=CSEL(13,I)/(CSEL(13,I)+CSIN(13,I))
      RATP =CSEL(14,I)/(CSEL(14,I)+CSIN(14,I))
      RATPB=CSEL(15,I)/(CSEL(15,I)+CSIN(15,I))
C     CSEL(11,I)=(6.*PAM           -3.*SAM   )  * RATKP
C     CSIN(11,I)=(6.*PAM           -3.*SAM   )  *(1.-RATKP)
C     CSEL(12,I)=(6.*PAM  +   AAM  -3.*SAM   )  * RATKM
C     CSIN(12,I)=(6.*PAM  +   AAM  -3.*SAM   )  *(1.-RATKM)
      CSEL(11,I)= CSEL(10,I)
      CSIN(11,I)= CSIN(10,I)
      CSEL(12,I)= CSEL(10,I)/4.+CSEL(12,I)/3.+5.*CSEL(13,I)/12.
      CSIN(12,I)= CSIN(10,I)/4.+CSIN(12,I)/3.+5.*CSIN(13,I)/12.
      PAMP= (CSEL(14,I)+CSIN(14,I))/9.
      AAMP= (CSEL(15,I)+CSIN(15,I)-CSEL(14,I)-CSIN(14,I))/5.
      SAMP= 2.*(CSEL(14,I)+CSIN(14,I))/9.-(CSEL(10,I)+CSIN(10,I))/3.
      IF(PLAB(I).LT.0.59) SAMP = 0.
      IF(PRFLAG) PRINT 1009,PLAB(I),PAM,AAM,SAM,PAMP,AAMP,SAMP
 1009 FORMAT(1H ,7F10.2)
      CSEL(17,I)=(9.*PAMP +4.*AAMP           )  *RATPB
      CSIN(17,I)=(9.*PAMP +4.*AAMP           )  *(1.-RATPB)
      CSEL(18,I)=(9.*PAMP          -3.*SAMP  )  *RATP
      CSIN(18,I)=(9.*PAMP          -3.*SAMP  )  *(1.-RATP)
      CSEL(19,I)=(9.*PAMP +2.*AAMP -3.*SAMP  )  *RATPB
      CSIN(19,I)=(9.*PAMP +2.*AAMP -3.*SAMP  )  *(1.-RATPB)
      CSEL(20,I)=(9.*PAMP          -3.*SAMP  )  *RATP
      CSIN(20,I)=(9.*PAMP          -3.*SAMP  )  *(1.-RATP)
      CSEL(22,I)=(9.*PAMP          -3.*SAMP  )  *RATP
      CSIN(22,I)=(9.*PAMP          -3.*SAMP  )  *(1.-RATP)
      CSEL(23,I)=(9.*PAMP +4.*AAMP -3.*SAMP  )  *RATPB
      CSIN(23,I)=(9.*PAMP +4.*AAMP -3.*SAMP  )  *(1.-RATPB)
      CSEL(25,I)=(9.*PAMP +2.*AAMP -3.*SAMP  )  *RATPB
      CSIN(25,I)=(9.*PAMP +2.*AAMP -3.*SAMP  )  *(1.-RATPB)
      CSEL(26,I)=(9.*PAMP          -6.*SAMP  )  *RATP
      CSIN(26,I)=(9.*PAMP          -6.*SAMP  )  *(1.-RATP)
      CSEL(27,I)=(9.*PAMP          -6.*SAMP  )  *RATP
      CSIN(27,I)=(9.*PAMP          -6.*SAMP  )  *(1.-RATP)
      CSEL(28,I)=(9.*PAMP +2.*AAMP -6.*SAMP  )  *RATPB
      CSIN(28,I)=(9.*PAMP +2.*AAMP -6.*SAMP  )  *(1.-RATPB)
      CSEL(29,I)=(9.*PAMP +   AAMP -6.*SAMP  )  *RATPB
      CSIN(29,I)=(9.*PAMP +   AAMP -6.*SAMP  )  *(1.-RATPB)
C --- SET CROSS SECTIONS FOR THE OMEGA AND ANTI-OMEGA TO THE VALUES ---
C --- AS FOR XI- AND ANTI XI- RESPECTIVELY ---
      CSEL(33,I)=CSEL(27,I)
      CSIN(33,I)=CSIN(27,I)
      CSEL(34,I)=CSEL(29,I)
      CSIN(34,I)=CSIN(29,I)
 41   CONTINUE
      DO 42 I=1,41
      DO 42 J=1,35
      IF(CSEL(J,I).LT.0.) CSEL(J,I)=0.
      IF(CSIN(J,I).LT.0.) CSIN(J,I)=0.
   42 CONTINUE
      GOTO 1
   60 DO 61 I=1,41
      CSPIEL(1,I)=A(1,I)-A(2,I)
      CSPIIN(1,I)=A(2,I)
      CSPIEL(2,I)=A(3,I)-A(4,I)
      CSPIIN(2,I)=A(4,I)
      CSPIEL(3,I)=A(5,I)-A(6,I)
   61 CSPIIN(3,I)=A(6,I)
      GOTO 1
   70 DO 71 I=1,41
      CSPNEL(1,I)=A(1,I)-A(2,I)
      CSPNIN(1,I)=A(2,I)
      CSPNEL(2,I)=A(3,I)-A(4,I)
      CSPNIN(2,I)=A(4,I)
      CSPNEL(3,I)=A(5,I)-A(6,I)
   71 CSPNIN(3,I)=A(6,I)
      GOTO 1
   80 IAT=IAT+3
      DO 81 I=1,17
      CNLWEL(IAT+1,I)=A(1,I)-A(2,I)
      CNLWIN(IAT+1,I)=A(2,I)
      CNLWEL(IAT+2,I)=A(3,I)-A(4,I)
      CNLWIN(IAT+2,I)=A(4,I)
      CNLWEL(IAT+3,I)=A(5,I)-A(6,I)
   81 CNLWIN(IAT+3,I)=A(6,I)
      GOTO 1
   90 I2=0
      I1=0
      DO 91 I=1,17
      I1=I2+1
      I2=I2+6
      IF(I2.GT.100) I2=100
      READ(LUNIN,1011) IA,(CSCAP(J),J=I1,I2)
      IF(PRFLAG)
     *PRINT 1011,IA,(CSCAP(J),J=I1,I2)
 1011 FORMAT(A10,6F10.2)
   91 CONTINUE
      GOTO 1
   95 READ(LUNIN,1001) IA
      IF(PRFLAG) PRINT 1002,IA
      DO 96 I=1,21
      READ(LUNIN,1012) EKFISS(I),(CSFISS(J,I),J=1,4)
      IF(PRFLAG) PRINT 1012,EKFISS(I),(CSFISS(J,I),J=1,4)
 1012 FORMAT(1X,F9.4,4F10.0)
   96 CONTINUE
C**
  100 IF(.NOT.PRFLAG) GOTO 200
      PRINT 1005
 1005 FORMAT(1H0,'LISTINGS OF ALL PARTICLE CROSS SECTIONS ON PROTONS MEA
     *SURED OR CALCULATED')
      DO 105 K=1,7
      K1=(K-1)*5+1
      K2=K*5
      PRINT 1006,K1,K2
 1006 FORMAT(1H0,'IPART1 - IPART2',2X,I3,'-',I3)
      DO 104 I=1,41
      PRINT 1007,PLAB(I),(CSEL(J,I),CSIN(J,I),J=K1,K2)
 1007 FORMAT(1H ,F10.2,2X,5(2X,2F10.2))
  104 CONTINUE
  105 CONTINUE
C
  200 CONTINUE
C
C --- CHECK FOR UN-PHYSICAL VALUES ---
      DO 9000 J=1,41
      IF (PLAB(J) .LT. 0.0) PLAB(J)=0.0
      DO 9001 I=1,35
      IF (CSEL(I,J) .LT. 0.0) CSEL(I,J)=0.0
      IF (CSIN(I,J) .LT. 0.0) CSIN(I,J)=0.0
      IF (I .GT. 3) GO TO 9001
      IF (CSPIEL(I,J) .LT. 0.0) CSPIEL(I,J)=0.0
      IF (CSPIIN(I,J) .LT. 0.0) CSPIIN(I,J)=0.0
      IF (CSPNEL(I,J) .LT. 0.0) CSPNEL(I,J)=0.0
      IF (CSPNIN(I,J) .LT. 0.0) CSPNIN(I,J)=0.0
 9001 CONTINUE
 9000 CONTINUE
C
      DO 9002 J=1,17
      IF (ELAB(J) .LT. 0.0) ELAB(J)=0.0
      DO 9003 I=1,15
      IF (CNLWEL(I,J) .LT. 0.0) CNLWEL(I,J)=0.0
      IF (CNLWIN(I,J) .LT. 0.0) CNLWIN(I,J)=0.0
 9003 CONTINUE
      IF (J .GT. 15) GO TO 9002
      IF (CNLWAT(J) .LT. 0.0) CNLWAT(J)=0.0
 9002 CONTINUE
C
      DO 9004 J=1,21
      IF (EKFISS(J) .LT. 0.0) EKFISS(J)=0.0
      DO 9005 I=1,4
      IF (CSFISS(I,J) .LT. 0.0) CSFISS(I,J)=0.0
 9005 CONTINUE
 9004 CONTINUE
C
      DO 9006 I=1,100
      IF (CSCAP(I) .LT. 0.0) CSCAP(I)=0.0
 9006 CONTINUE
C
C --- WRITE OUT THE CROSS SECTIONS IN THE FORM OF DATA STMTS. ---
      WRITE(LUNOUT,2001)
 2001 FORMAT('+KEEP,/CSDAT.'/
     $ 'C --- CROSS-SECTION DATA BY "PCSDAT" 01-FEB-1989 ---'/'C')
C
      WRITE(LUNOUT,2010) (PLAB(I),I=1,41)
 2010 FORMAT(6X,'DATA PLAB /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
C
      DO 220 I=1,35
      WRITE(LUNOUT,2020) I,(CSEL(I,J),J=1,41)
 2020 FORMAT(6X,'DATA (CSEL(',I2,',J),J=1,41) /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 220  CONTINUE
C
      DO 230 I=1,35
      WRITE(LUNOUT,2030) I,(CSIN(I,J),J=1,41)
 2030 FORMAT(6X,'DATA (CSIN(',I2,',J),J=1,41) /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 230  CONTINUE
C
      DO 240 I=1,3
      WRITE(LUNOUT,2040) I,(CSPIEL(I,J),J=1,41)
 2040 FORMAT(6X,'DATA (CSPIEL(',I2,',J),J=1,41) /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 240  CONTINUE
C
      DO 250 I=1,3
      WRITE(LUNOUT,2050) I,(CSPIIN(I,J),J=1,41)
 2050 FORMAT(6X,'DATA (CSPIIN(',I2,',J),J=1,41) /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 250  CONTINUE
C
      DO 260 I=1,3
      WRITE(LUNOUT,2060) I,(CSPNEL(I,J),J=1,41)
 2060 FORMAT(6X,'DATA (CSPNEL(',I2,',J),J=1,41) /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 260  CONTINUE
C
      DO 270 I=1,3
      WRITE(LUNOUT,2070) I,(CSPNIN(I,J),J=1,41)
 2070 FORMAT(6X,'DATA (CSPNIN(',I2,',J),J=1,41) /'/
     $ 8(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 270  CONTINUE
C
      WRITE(LUNOUT,2080) (ELAB(I),I=1,17)
 2080 FORMAT(6X,'DATA ELAB /'/
     $ 3(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,',',G12.5,'/')
C
      WRITE(LUNOUT,2090) (CNLWAT(I),I=1,15)
 2090 FORMAT(6X,'DATA CNLWAT /'/
     $ 2(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',4(G12.5,','),G12.5,'/')
C
      DO 280 I=1,15
      WRITE(LUNOUT,2100) I,(CNLWEL(I,J),J=1,17)
 2100 FORMAT(6X,'DATA (CNLWEL(',I2,',J),J=1,17) /'/
     $ 3(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,',',G12.5,'/')
 280  CONTINUE
C
      DO 290 I=1,15
      WRITE(LUNOUT,2110) I,(CNLWIN(I,J),J=1,17)
 2110 FORMAT(6X,'DATA (CNLWIN(',I2,',J),J=1,17) /'/
     $ 3(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,',',G12.5,'/')
 290  CONTINUE
C
C --- WRITE CSCAP IN TWO PARTS BECAUSE OF CONTINUATION CARD LIMIT ---
      WRITE(LUNOUT,2120) (CSCAP(I),I=1,100)
 2120 FORMAT(6X,'DATA (CSCAP(J),J=1,50) /'/
     $ 9(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',4(G12.5,','),G12.5,'/'/
     $       6X,'DATA (CSCAP(J),J=51,100) /'/
     $ 9(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',4(G12.5,','),G12.5,'/')
C
      WRITE(LUNOUT,2130) (EKFISS(I),I=1,21)
 2130 FORMAT(6X,'DATA EKFISS /'/
     $ 4(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
C
      DO 300 I=1,4
      WRITE(LUNOUT,2140) I,(CSFISS(I,J),J=1,21)
 2140 FORMAT(6X,'DATA (CSFISS(',I2,',J),J=1,21) /'/
     $ 4(5X,'$ ',5(G12.5,',')/),
     $   5X,'$ ',G12.5,'/')
 300  CONTINUE
C
      WRITE(LUNOUT,2999)
 2999 FORMAT('C --- END OF CROSS SECTION DATA STATEMENTS ---'/'C')
C
      RETURN
      END

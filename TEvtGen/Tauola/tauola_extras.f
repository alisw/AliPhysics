
      SUBROUTINE CHOICE(MNUM,RR,ICHAN,PROB1,PROB2,PROB3,
     $            AMRX,GAMRX,AMRA,GAMRA,AMRB,GAMRB)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      AMROP=1.1
      GAMROP=0.36
      AMOM=.782
      GAMOM=0.0084
C     XXXXA CORRESPOND TO S2 CHANNEL !
      IF(MNUM.EQ.0) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =AMA1
       GAMRX=GAMA1
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMRO
       GAMRB=GAMRO
      ELSEIF(MNUM.EQ.1) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.57
       GAMRX=0.9
       AMRB =AMKST
       GAMRB=GAMKST
       AMRA =AMRO
       GAMRA=GAMRO
      ELSEIF(MNUM.EQ.2) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.57
       GAMRX=0.9
       AMRB =AMKST
       GAMRB=GAMKST
       AMRA =AMRO
       GAMRA=GAMRO
      ELSEIF(MNUM.EQ.3) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMKST
       GAMRA=GAMKST
       AMRB =AMKST
       GAMRB=GAMKST
      ELSEIF(MNUM.EQ.4) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMKST
       GAMRA=GAMKST
       AMRB =AMKST
       GAMRB=GAMKST
      ELSEIF(MNUM.EQ.5) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMKST
       GAMRA=GAMKST
       AMRB =AMRO
       GAMRB=GAMRO
      ELSEIF(MNUM.EQ.6) THEN
       PROB1=0.4
       PROB2=0.4
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMKST
       GAMRB=GAMKST
      ELSEIF(MNUM.EQ.7) THEN
       PROB1=0.0
       PROB2=1.0
       AMRX =1.27
       GAMRX=0.9
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMRO
       GAMRB=GAMRO
      ELSEIF(MNUM.EQ.8) THEN
       PROB1=0.0
       PROB2=1.0
       AMRX =AMROP
       GAMRX=GAMROP
       AMRB =AMOM
       GAMRB=GAMOM
       AMRA =AMRO
       GAMRA=GAMRO
      ELSEIF(MNUM.EQ.101) THEN
       PROB1=.35
       PROB2=.35
       AMRX =1.2
       GAMRX=.46
       AMRB =AMOM
       GAMRB=GAMOM
       AMRA =AMOM
       GAMRA=GAMOM
      ELSEIF(MNUM.EQ.102) THEN
       PROB1=0.0
       PROB2=0.0
       AMRX =1.4
       GAMRX=.6
       AMRB =AMOM
       GAMRB=GAMOM
       AMRA =AMOM
       GAMRA=GAMOM
      ELSE
       PROB1=0.0
       PROB2=0.0
       AMRX =AMA1
       GAMRX=GAMA1
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMRO
       GAMRB=GAMRO
      ENDIF
C
      IF    (RR.LE.PROB1) THEN
       ICHAN=1
      ELSEIF(RR.LE.(PROB1+PROB2)) THEN
       ICHAN=2
        AX   =AMRA
        GX   =GAMRA
        AMRA =AMRB
        GAMRA=GAMRB
        AMRB =AX
        GAMRB=GX
        PX   =PROB1
        PROB1=PROB2
        PROB2=PX
      ELSE
       ICHAN=3
      ENDIF
C
      PROB3=1.0-PROB1-PROB2
      END
      SUBROUTINE INITDK
* ----------------------------------------------------------------------
*     INITIALISATION OF TAU DECAY PARAMETERS  and routines
*
*     called by : KORALZ
* ----------------------------------------------------------------------

      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS






      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31

      CHARACTER OLDNAMES(7)*31
      CHARACTER*80 bxINIT
      PARAMETER (
     $  bxINIT ='(1x,1h*,g17.8,            16x, a31,a4,a4, 1x,1h*)'
     $ )
      REAL*4 PI
*
*
* LIST OF BRANCHING RATIOS
CAM normalised to e nu nutau channel
CAM                  enu   munu   pinu  rhonu   A1nu   Knu    K*nu   pi
CAM   DATA JLIST  /    1,     2,     3,     4,     5,     6,     7,
*AM   DATA GAMPRT /1.000,0.9730,0.6054,1.2432,0.8432,0.0432,O.O811,0.616
*AM
*AM  multipion decays
*
*    conventions of particles names
*                 K-,P-,K+,  K0,P-,KB,  K-,P0,K0
*                  3, 1,-3  , 4, 1,-4  , 3, 2, 4  ,
*                 P0,P0,K-,  K-,P-,P+,  P-,KB,P0
*                  2, 2, 3  , 3, 1,-1  , 1,-4, 2  ,
*                 ET,P-,P0   P-,P0,GM
*                  9, 1, 2  , 1, 2, 8
*
C
      DIMENSION NOPIK(6,NMODE),NPIK(NMODE)
*AM   outgoing multiplicity and flavors of multi-pion /multi-K modes    
      DATA   NPIK  /                4,                    4,  
     1                              5,                    5,
     2                              6,                    6,
     3                              3,                    3,            
     4                              3,                    3,            
     5                              3,                    3,            
     6                              3,                    3,  
     7                              2                         /         
      DATA  NOPIK / -1,-1, 1, 2, 0, 0,     2, 2, 2,-1, 0, 0,  
     1              -1,-1, 1, 2, 2, 0,    -1,-1,-1, 1, 1, 0,  
     2              -1,-1,-1, 1, 1, 2,    -1,-1, 1, 2, 2, 2, 
     3              -3,-1, 3, 0, 0, 0,    -4,-1, 4, 0, 0, 0,  
     4              -3, 2,-4, 0, 0, 0,     2, 2,-3, 0, 0, 0,  
     5              -3,-1, 1, 0, 0, 0,    -1, 4, 2, 0, 0, 0,  
     6               9,-1, 2, 0, 0, 0,    -1, 2, 8, 0, 0, 0,
C AJWMOD fix sign bug, 2/22/99
     7              -3,-4, 0, 0, 0, 0                         /
* LIST OF BRANCHING RATIOS
      NCHAN = NMODE + 7
      DO 1 I = 1,30
      IF (I.LE.NCHAN) THEN
        JLIST(I) = I
        IF(I.EQ. 1) GAMPRT(I) =0.1800 
        IF(I.EQ. 2) GAMPRT(I) =0.1751 
        IF(I.EQ. 3) GAMPRT(I) =0.1110 
        IF(I.EQ. 4) GAMPRT(I) =0.2515 
        IF(I.EQ. 5) GAMPRT(I) =0.1790 
        IF(I.EQ. 6) GAMPRT(I) =0.0071 
        IF(I.EQ. 7) GAMPRT(I) =0.0134
        IF(I.EQ. 8) GAMPRT(I) =0.0450
        IF(I.EQ. 9) GAMPRT(I) =0.0100
        IF(I.EQ.10) GAMPRT(I) =0.0009
        IF(I.EQ.11) GAMPRT(I) =0.0004 
        IF(I.EQ.12) GAMPRT(I) =0.0003 
        IF(I.EQ.13) GAMPRT(I) =0.0005 
        IF(I.EQ.14) GAMPRT(I) =0.0015 
        IF(I.EQ.15) GAMPRT(I) =0.0015 
        IF(I.EQ.16) GAMPRT(I) =0.0015 
        IF(I.EQ.17) GAMPRT(I) =0.0005
        IF(I.EQ.18) GAMPRT(I) =0.0050
        IF(I.EQ.19) GAMPRT(I) =0.0055
        IF(I.EQ.20) GAMPRT(I) =0.0017 
        IF(I.EQ.21) GAMPRT(I) =0.0013 
        IF(I.EQ.22) GAMPRT(I) =0.0010 
        IF(I.EQ. 1) OLDNAMES(I)='  TAU-  -->   E-               '
        IF(I.EQ. 2) OLDNAMES(I)='  TAU-  -->  MU-               '
        IF(I.EQ. 3) OLDNAMES(I)='  TAU-  -->  PI-               '
        IF(I.EQ. 4) OLDNAMES(I)='  TAU-  -->  PI-, PI0          '
        IF(I.EQ. 5) OLDNAMES(I)='  TAU-  -->  A1- (two subch)   '
        IF(I.EQ. 6) OLDNAMES(I)='  TAU-  -->   K-               '
        IF(I.EQ. 7) OLDNAMES(I)='  TAU-  -->  K*- (two subch)   '
        IF(I.EQ. 8) NAMES(I-7)='  TAU-  --> 2PI-,  PI0,  PI+   '
        IF(I.EQ. 9) NAMES(I-7)='  TAU-  --> 3PI0,        PI-   '
        IF(I.EQ.10) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 2PI0   '
        IF(I.EQ.11) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,        '
        IF(I.EQ.12) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,  PI0   '
        IF(I.EQ.13) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 3PI0   '
        IF(I.EQ.14) NAMES(I-7)='  TAU-  -->  K-, PI-,  K+      '
        IF(I.EQ.15) NAMES(I-7)='  TAU-  -->  K0, PI-, K0B      '
        IF(I.EQ.16) NAMES(I-7)='  TAU-  -->  K-,  K0, PI0      '
        IF(I.EQ.17) NAMES(I-7)='  TAU-  --> PI0  PI0   K-      '
        IF(I.EQ.18) NAMES(I-7)='  TAU-  -->  K-  PI-  PI+      '
        IF(I.EQ.19) NAMES(I-7)='  TAU-  --> PI-  K0B  PI0      '
        IF(I.EQ.20) NAMES(I-7)='  TAU-  --> ETA  PI-  PI0      '
        IF(I.EQ.21) NAMES(I-7)='  TAU-  --> PI-  PI0  GAM      '
        IF(I.EQ.22) NAMES(I-7)='  TAU-  -->  K-  K0            '
      ELSE
        JLIST(I) = 0
        GAMPRT(I) = 0.
      ENDIF
   1  CONTINUE
      DO I=1,NMODE
        MULPIK(I)=NPIK(I)
        DO J=1,MULPIK(I)
         IDFFIN(J,I)=NOPIK(J,I)
        ENDDO
      ENDDO
*
*
* --- COEFFICIENTS TO FIX RATIO OF:
* --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
* --- PROBABILITY OF K0 TO BE KS
* --- PROBABILITY OF K0B TO BE KS
* --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
* --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
* --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
* --- NEGLECTS MASS-PHASE SPACE EFFECTS
      BRA1=0.5
      BRK0=0.5
      BRK0B=0.5
      BRKS=0.6667
*

      GFERMI = 1.16637E-5
      CCABIB = 0.975
      GV     = 1.0
      GA     =-1.0



* ZW 13.04.89 HERE WAS AN ERROR
      SCABIB = SQRT(1.-CCABIB**2)
      PI =4.*ATAN(1.)
      GAMEL  = GFERMI**2*AMTAU**5/(192*PI**3)
*
*      CALL DEXAY(-1,pol1)
*
      RETURN
      END
      FUNCTION DCDMAS(IDENT)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      IF      (IDENT.EQ. 1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ.-1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ. 2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ.-2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ. 3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ.-3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ. 4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ.-4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ. 8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ.-8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ. 9) THEN
        APKMAS=0.5488
      ELSEIF  (IDENT.EQ.-9) THEN
        APKMAS=0.5488
      ELSE
        PRINT *, 'STOP IN APKMAS, WRONG IDENT=',IDENT
        STOP
      ENDIF
      DCDMAS=APKMAS
      END
      FUNCTION LUNPIK(ID,ISGN)
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      REAL*4 XIO(1)
      IDENT=ID*ISGN
      IF      (IDENT.EQ. 1) THEN
        IPKDEF=-211
      ELSEIF  (IDENT.EQ.-1) THEN
        IPKDEF= 211
      ELSEIF  (IDENT.EQ. 2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ.-2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ. 3) THEN
        IPKDEF=-321
      ELSEIF  (IDENT.EQ.-3) THEN
        IPKDEF= 321
      ELSEIF  (IDENT.EQ. 4) THEN
*
* K0 --> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO(1).GT.BRK0) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ.-4) THEN
*
* K0B--> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO(1).GT.BRK0B) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ. 8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ.-8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ. 9) THEN
        IPKDEF= 221
      ELSEIF  (IDENT.EQ.-9) THEN
        IPKDEF= 221
      ELSE
        PRINT *, 'STOP IN IPKDEF, WRONG IDENT=',IDENT
        STOP
      ENDIF
      LUNPIK=IPKDEF
      END



      SUBROUTINE TAURDF(KTO)
C THIS ROUTINE CAN BE CALLED BEFORE ANY TAU+ OR TAU- EVENT IS GENERATED
C IT CAN BE USED TO GENERATE TAU+ AND TAU- SAMPLES OF DIFFERENT
C CONTENTS
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN

C Subroutine TAURDF is disabled
      RETURN


      IF (KTO.EQ.1) THEN
C     ==================
C AJWMOD: Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
      BRA1 = PKORB(4,1)
      BRKS = PKORB(4,3)
      BRK0  = PKORB(4,5)
      BRK0B  = PKORB(4,6)
      ELSE
C     ====
C AJWMOD: Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
      BRA1 = PKORB(4,2)
      BRKS = PKORB(4,4)
      BRK0  = PKORB(4,5)
      BRK0B  = PKORB(4,6)
      ENDIF
C     =====
      END

      SUBROUTINE INIPHY(XK00)
* ----------------------------------------------------------------------
*     INITIALISATION OF PARAMETERS
*     USED IN QED and/or GSW ROUTINES
* ----------------------------------------------------------------------
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      REAL*8 PI8,XK00
*
      PI8    = 4.D0*DATAN(1.D0)
      ALFINV = 137.03604D0
      ALFPI  = 1D0/(ALFINV*PI8)
      XK0=XK00
      END

      SUBROUTINE INIMAS
C ----------------------------------------------------------------------
C     INITIALISATION OF MASSES
C
C     called by : KORALZ
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
C IN-COMING / OUT-GOING  FERMION MASSES
      AMTAU  = 1.7842
C --- let us update tau mass ...
      AMTAU  = 1.777
      AMNUTA = 0.010
      AMEL   = 0.0005111
      AMNUE  = 0.0
      AMMU   = 0.105659 
      AMNUMU = 0.0
*
* MASSES USED IN TAU DECAYS
      AMPIZ  = 0.134964
      AMPI   = 0.139568
      AMRO   = 0.773
      GAMRO  = 0.145
*C    GAMRO  = 0.666
      AMA1   = 1.251
      GAMA1  = 0.599
      AMK    = 0.493667
      AMKZ   = 0.49772
      AMKST  = 0.8921
      GAMKST = 0.0513
C
C
C IN-COMING / OUT-GOING  FERMION MASSES
!!      AMNUTA = PKORB(1,2)
!!      AMNUE  = PKORB(1,4)
!!      AMNUMU = PKORB(1,6)
C
C MASSES USED IN TAU DECAYS  Cleo settings
!!      AMPIZ  = PKORB(1,7)
!!      AMPI   = PKORB(1,8)
!!      AMRO   = PKORB(1,9)
!!      GAMRO  = PKORB(2,9)
      AMA1   = 1.275   !! PKORB(1,10)
      GAMA1  = 0.615   !! PKORB(2,10)
!!      AMK    = PKORB(1,11)
!!      AMKZ   = PKORB(1,12)
!!      AMKST  = PKORB(1,13)
!!      GAMKST = PKORB(2,13)
C

      RETURN
      END
      SUBROUTINE ANGULU(PD1,PD2,Q1,Q2,COSTHE)
      REAL*8 PD1(4),PD2(4),Q1(4),Q2(4),COSTHE,P(4),QQ(4),QT(4)
C take effective beam which is less massive, it should be irrelevant
C but in case HEPEVT is particulary dirty may help.
C this routine calculate reduced system transver and cosine of scattering 
C angle.

      XM1=ABS(PD1(4)**2-PD1(3)**2-PD1(2)**2-PD1(1)**2)
      XM2=ABS(PD2(4)**2-PD2(3)**2-PD2(2)**2-PD2(1)**2)
      IF (XM1.LT.XM2) THEN
        SIGN=1D0
        DO K=1,4
          P(K)=PD1(K)
        ENDDO
      ELSE
        SIGN=-1D0
        DO K=1,4
          P(K)=PD2(K)
        ENDDO
      ENDIF
C calculate space like part of P (in Z restframe)
      DO K=1,4
       QQ(K)=Q1(k)+Q2(K)
       QT(K)=Q1(K)-Q2(K)
      ENDDO

       XMQQ=SQRT(QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2)

       QTXQQ=QT(4)*QQ(4)-QT(3)*QQ(3)-QT(2)*QQ(2)-QT(1)*QQ(1)
      DO K=1,4
       QT(K)=QT(K)-QQ(K)*QTXQQ/XMQQ**2
      ENDDO

       PXQQ=P(4)*QQ(4)-P(3)*QQ(3)-P(2)*QQ(2)-P(1)*QQ(1)
      DO K=1,4
       P(K)=P(K)-QQ(K)*PXQQ/XMQQ**2
      ENDDO
C calculate costhe
       PXP  =SQRT(p(1)**2+p(2)**2+p(3)**2-p(4)**2)
       QTXQT=SQRT(QT(3)**2+QT(2)**2+QT(1)**2-QT(4)**2)
       PXQT =P(3)*QT(3)+P(2)*QT(2)+P(1)*QT(1)-P(4)*QT(4)
       COSTHE=PXQT/PXP/QTXQT
       COSTHE=COSTHE*SIGN
      END

      FUNCTION PLZAP0(IDE,IDF,SVAR,COSTH0)
C this function calculates probability for the helicity +1 +1 configuration
C of taus for given Z/gamma transfer and COSTH0 cosine of scattering angle
      REAL*8 PLZAP0,SVAR,COSTHE,COSTH0,T_BORN

      COSTHE=COSTH0
C >>>>>      IF (IDE*IDF.LT.0) COSTHE=-COSTH0 ! this is probably not needed ID
C >>>>>      of first beam is used by T_GIVIZ0 including sign

      IF (IDF.GT.0) THEN
        CALL INITWK(IDE,IDF,SVAR)
      ELSE
        CALL INITWK(-IDE,-IDF,SVAR)
      ENDIF
      PLZAP0=T_BORN(0,SVAR,COSTHE,1D0,1D0)
     $  /(T_BORN(0,SVAR,COSTHE,1D0,1D0)+T_BORN(0,SVAR,COSTHE,-1D0,-1D0))


C      write(*,*) 'svar=  ',  svar
C      write(*,*) 'COSTHE=',  COSTHE
C      write(*,*) ide,'  ',idf
C      write(*,*) 'PLZAP0=', PLZAP0
C      COSTHE=0.999
C      write(*,*) 'TBORN+=', T_BORN(0,SVAR,COSTHE,1D0,1D0)

C      COSTHE=-0.999
C      write(*,*) 'TBORN-=', T_BORN(0,SVAR,COSTHE,-1D0,-1D0)
C 100  format (A,E8.4)

!      PLZAP0=0.5
      END
      FUNCTION T_BORN(MODE,SVAR,COSTHE,TA,TB)
C ----------------------------------------------------------------------
C THIS ROUTINE PROVIDES BORN CROSS SECTION. IT HAS THE SAME         
C STRUCTURE AS FUNTIS AND FUNTIH, THUS CAN BE USED AS SIMPLER       
C EXAMPLE OF THE METHOD APPLIED THERE                               
C INPUT PARAMETERS ARE: SVAR    -- transfer
C                       COSTHE  -- cosine of angle between tau+ and 1st beam
C                       TA,TB   -- helicity states of tau+ tau-
C
C     called by : BORNY, BORAS, BORNV, WAGA, WEIGHT
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / T_BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF
      REAL*8              ENE ,AMIN,AMFIN
      COMMON / T_GAUSPM /SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI   ,XUPZI   ,XUPGF   ,XUPZF
     &                  ,NDIAG0,NDIAGA,KEYA,KEYZ
     &                  ,ITCE,JTCE,ITCF,JTCF,KOLOR
      REAL*8             SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI(2),XUPZI(2),XUPGF(2),XUPZF(2)
      REAL*8            SEPS1,SEPS2
C=====================================================================
      COMMON / T_GSWPRM /SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      REAL*8             SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
C     SWSQ        = sin2 (theta Weinberg)
C     AMW,AMZ     = W & Z boson masses respectively
C     AMH         = the Higgs mass
C     AMTOP       = the top mass
C     GAMMZ       = Z0 width
      COMPLEX*16 ABORN(2,2),APHOT(2,2),AZETT(2,2)
      COMPLEX*16 XUPZFP(2),XUPZIP(2)
      COMPLEX*16 ABORNM(2,2),APHOTM(2,2),AZETTM(2,2)
      COMPLEX*16 PROPA,PROPZ
      COMPLEX*16 XR,XI
      COMPLEX*16 XUPF,XUPI
      COMPLEX*16 XTHING
      DATA XI/(0.D0,1.D0)/,XR/(1.D0,0.D0)/
      DATA MODE0 /-5/
      DATA IDE0 /-55/
      DATA SVAR0,COST0 /-5.D0,-6.D0/
      DATA PI /3.141592653589793238462643D0/
      DATA SEPS1,SEPS2 /0D0,0D0/
C
C MEMORIZATION =========================================================
      IF ( MODE.NE.MODE0.OR.SVAR.NE.SVAR0.OR.COSTHE.NE.COST0
     $    .OR.IDE0.NE.IDE)THEN
C
        KEYGSW=1
C ** PROPAGATORS
        IDE0=IDE
        MODE0=MODE
        SVAR0=SVAR
        COST0=COSTHE
        SINTHE=SQRT(1.D0-COSTHE**2)
        BETA=SQRT(MAX(0D0,1D0-4D0*AMFIN**2/SVAR))
C I MULTIPLY AXIAL COUPLING BY BETA FACTOR.
        XUPZFP(1)=0.5D0*(XUPZF(1)+XUPZF(2))+0.5*BETA*(XUPZF(1)-XUPZF(2))
        XUPZFP(2)=0.5D0*(XUPZF(1)+XUPZF(2))-0.5*BETA*(XUPZF(1)-XUPZF(2))
        XUPZIP(1)=0.5D0*(XUPZI(1)+XUPZI(2))+0.5*(XUPZI(1)-XUPZI(2))
        XUPZIP(2)=0.5D0*(XUPZI(1)+XUPZI(2))-0.5*(XUPZI(1)-XUPZI(2))
C FINAL STATE VECTOR COUPLING
        XUPF     =0.5D0*(XUPZF(1)+XUPZF(2))
        XUPI     =0.5D0*(XUPZI(1)+XUPZI(2))
        XTHING   =0D0

        PROPA =1D0/SVAR
        PROPZ =1D0/DCMPLX(SVAR-AMZ**2,SVAR/AMZ*GAMMZ)
        IF (KEYGSW.EQ.0) PROPZ=0.D0
        DO 50 I=1,2
         DO 50 J=1,2
          REGULA= (3-2*I)*(3-2*J) + COSTHE
          REGULM=-(3-2*I)*(3-2*J) * SINTHE *2.D0*AMFIN/SQRT(SVAR)
          APHOT(I,J)=PROPA*(XUPGI(I)*XUPGF(J)*REGULA)
          AZETT(I,J)=PROPZ*(XUPZIP(I)*XUPZFP(J)+XTHING)*REGULA
          ABORN(I,J)=APHOT(I,J)+AZETT(I,J)
          APHOTM(I,J)=PROPA*DCMPLX(0D0,1D0)*XUPGI(I)*XUPGF(J)*REGULM
          AZETTM(I,J)=PROPZ*DCMPLX(0D0,1D0)*(XUPZIP(I)*XUPF+XTHING)*REGULM
          ABORNM(I,J)=APHOTM(I,J)+AZETTM(I,J)
   50   CONTINUE
      ENDIF
C
C******************
C* IN CALCULATING CROSS SECTION ONLY DIAGONAL ELEMENTS
C* OF THE SPIN DENSITY MATRICES ENTER (LONGITUD. POL. ONLY.)
C* HELICITY CONSERVATION EXPLICITLY OBEYED
      POLAR1=  (SEPS1)
      POLAR2= (-SEPS2)
      BORN=0D0
      DO 150 I=1,2
       HELIC= 3-2*I
       DO 150 J=1,2
        HELIT=3-2*J
        FACTOR=KOLOR*(1D0+HELIC*POLAR1)*(1D0-HELIC*POLAR2)/4D0
        FACTOM=FACTOR*(1+HELIT*TA)*(1-HELIT*TB)
        FACTOR=FACTOR*(1+HELIT*TA)*(1+HELIT*TB)

        BORN=BORN+CDABS(ABORN(I,J))**2*FACTOR
C      MASS TERM IN BORN
        IF (MODE.GE.1) THEN
         BORN=BORN+CDABS(ABORNM(I,J))**2*FACTOM
        ENDIF

  150 CONTINUE
C************
      FUNT=BORN
      IF(FUNT.LT.0.D0)  FUNT=BORN

C
      IF (SVAR.GT.4D0*AMFIN**2) THEN
C PHASE SPACE THRESHOLD FACTOR
        THRESH=SQRT(1-4D0*AMFIN**2/SVAR)
        T_BORN= FUNT*SVAR**2*THRESH
      ELSE
        THRESH=0.D0
        T_BORN=0.D0
      ENDIF
C ZW HERE WAS AN ERROR 19. 05. 1989
!      write(*,*) 'KKKK ',PROPA,PROPZ,XUPGI,XUPGF,XUPZI,XUPZF
!      write(*,*) 'KKKK X',svar,costhe,TA,TB,T_BORN
      END

      SUBROUTINE INITWK(IDEX,IDFX,SVAR)
! initialization routine coupling masses etc.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / T_BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF
      REAL*8              ENE ,AMIN,AMFIN
      COMMON / T_GAUSPM /SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI   ,XUPZI   ,XUPGF   ,XUPZF
     &                  ,NDIAG0,NDIAGA,KEYA,KEYZ
     &                  ,ITCE,JTCE,ITCF,JTCF,KOLOR
      REAL*8             SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI(2),XUPZI(2),XUPGF(2),XUPZF(2)
      COMMON / T_GSWPRM /SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      REAL*8             SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
C     SWSQ        = sin2 (theta Weinberg)
C     AMW,AMZ     = W & Z boson masses respectively
C     AMH         = the Higgs mass
C     AMTOP       = the top mass
C     GAMMZ       = Z0 width
C
      ENE=SQRT(SVAR)/2
      AMIN=0.511D-3
      SWSQ=0.23147
      AMZ=91.1882
      GAMMZ=2.4952
      IF     (IDFX.EQ. 15) then       
        IDF=2  ! denotes tau +2 tau-
        AMFIN=1.77703 !this mass is irrelevant if small, used in ME only
      ELSEIF (IDFX.EQ.-15) then
        IDF=-2  ! denotes tau -2 tau-
        AMFIN=1.77703 !this mass is irrelevant if small, used in ME only
      ELSE
        WRITE(*,*) 'INITWK: WRONG IDFX'
        STOP
      ENDIF

      IF     (IDEX.EQ. 11) then      !electron
        IDE= 2
        AMIN=0.511D-3
      ELSEIF (IDEX.EQ.-11) then      !positron
        IDE=-2
        AMIN=0.511D-3
      ELSEIF (IDEX.EQ. 13) then      !mu+
        IDE= 2
        AMIN=0.105659
      ELSEIF (IDEX.EQ.-13) then      !mu-
        IDE=-2
        AMIN=0.105659
      ELSEIF (IDEX.EQ.  1) then      !d
        IDE= 4
        AMIN=0.05
      ELSEIF (IDEX.EQ.- 1) then      !d~
        IDE=-4
        AMIN=0.05
      ELSEIF (IDEX.EQ.  2) then      !u
        IDE= 3
        AMIN=0.02
      ELSEIF (IDEX.EQ.- 2) then      !u~
        IDE=-3
        AMIN=0.02
      ELSEIF (IDEX.EQ.  3) then      !s
        IDE= 4
        AMIN=0.3
      ELSEIF (IDEX.EQ.- 3) then      !s~
        IDE=-4
        AMIN=0.3
      ELSEIF (IDEX.EQ.  4) then      !c
        IDE= 3
        AMIN=1.3
      ELSEIF (IDEX.EQ.- 4) then      !c~
        IDE=-3
        AMIN=1.3
      ELSEIF (IDEX.EQ.  5) then      !b
        IDE= 4
        AMIN=4.5
      ELSEIF (IDEX.EQ.- 5) then      !b~
        IDE=-4
        AMIN=4.5
      ELSEIF (IDEX.EQ.  12) then     !nu_e
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 12) then     !nu_e~
        IDE=-1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.  14) then     !nu_mu
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 14) then     !nu_mu~
        IDE=-1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.  16) then     !nu_tau
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 16) then     !nu_tau~
        IDE=-1
        AMIN=0.1D-3

      ELSE
        WRITE(*,*) 'INITWK: WRONG IDEX'
        STOP
      ENDIF

C ----------------------------------------------------------------------
C
C     INITIALISATION OF COUPLING CONSTANTS AND FERMION-GAMMA / Z0 VERTEX
C
C     called by : KORALZ
C ----------------------------------------------------------------------
      ITCE=IDE/IABS(IDE)
      JTCE=(1-ITCE)/2
      ITCF=IDF/IABS(IDF)
      JTCF=(1-ITCF)/2
      CALL T_GIVIZO( IDE, 1,AIZOR,QE,KDUMM)
      CALL T_GIVIZO( IDE,-1,AIZOL,QE,KDUMM)
      XUPGI(1)=QE
      XUPGI(2)=QE
      T3E    = AIZOL+AIZOR
      XUPZI(1)=(AIZOR-QE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      XUPZI(2)=(AIZOL-QE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      CALL T_GIVIZO( IDF, 1,AIZOR,QF,KOLOR)
      CALL T_GIVIZO( IDF,-1,AIZOL,QF,KOLOR)
      XUPGF(1)=QF
      XUPGF(2)=QF
      T3F    =  AIZOL+AIZOR
      XUPZF(1)=(AIZOR-QF*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      XUPZF(2)=(AIZOL-QF*SWSQ)/SQRT(SWSQ*(1-SWSQ))
C
      NDIAG0=2
      NDIAGA=11
      KEYA  = 1
      KEYZ  = 1
C
C
      RETURN
      END

      SUBROUTINE T_GIVIZO(IDFERM,IHELIC,SIZO3,CHARGE,KOLOR)
C ----------------------------------------------------------------------
C PROVIDES ELECTRIC CHARGE AND WEAK IZOSPIN OF A FAMILY FERMION
C IDFERM=1,2,3,4 DENOTES NEUTRINO, LEPTON, UP AND DOWN QUARK
C NEGATIVE IDFERM=-1,-2,-3,-4, DENOTES ANTIPARTICLE
C IHELIC=+1,-1 DENOTES RIGHT AND LEFT HANDEDNES ( CHIRALITY)
C SIZO3 IS THIRD PROJECTION OF WEAK IZOSPIN (PLUS MINUS HALF)
C AND CHARGE IS ELECTRIC CHARGE IN UNITS OF ELECTRON CHARGE
C KOLOR IS A QCD COLOUR, 1 FOR LEPTON, 3 FOR QUARKS
C
C     called by : EVENTE, EVENTM, FUNTIH, .....
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
C
      IF(IDFERM.EQ.0.OR.IABS(IDFERM).GT.4) GOTO 901
      IF(IABS(IHELIC).NE.1)                GOTO 901
      IH  =IHELIC
      IDTYPE =IABS(IDFERM)
      IC  =IDFERM/IDTYPE
      LEPQUA=INT(IDTYPE*0.4999999D0)
      IUPDOW=IDTYPE-2*LEPQUA-1
      CHARGE  =(-IUPDOW+2D0/3D0*LEPQUA)*IC
      SIZO3   =0.25D0*(IC-IH)*(1-2*IUPDOW)
      KOLOR=1+2*LEPQUA
C** NOTE THAT CONVENTIONALY Z0 COUPLING IS
C** XOUPZ=(SIZO3-CHARGE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
      RETURN
 901  PRINT *,' STOP IN GIVIZO: WRONG PARAMS.'
      STOP
      END
      SUBROUTINE PHYFIX(NSTOP,NSTART)
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5) 
      SAVE /LUJETS/ 
C NSTOP NSTART : when PHYTIA history ends and event starts.
      NSTOP=0
      NSTART=1
      DO I=1, N
       IF(K(I,1).NE.21) THEN
           NSTOP = I-1
           NSTART= I
           GOTO 500
       ENDIF
      ENDDO
 500  CONTINUE
      END
 

      SUBROUTINE TAUPI0(PI0,K)
C no initialization required. Must be called once after every:
C   1)    CALL DEKAY(1+10,...)
C   2)    CALL DEKAY(2+10,...)
C   3)    CALL DEXAY(1,...)
C   4)    CALL DEXAY(2,...)
C subroutine to decay originating from TAUOLA's taus: 
C 1) etas (with CALL TAUETA(JAK))
C 2) later pi0's from taus.
C 3) extensions to other applications possible. 
C this routine belongs to >tauola universal interface<, but uses 
C routines from >tauola< utilities as well.  25.08.2005      
C this is the hepevt class in old style. No d_h_ class pre-name
 
C position of taus, must be defined by host program:
      COMMON /TAUPOS/ NP1,NP2
c
      REAL  PHOT1(4),PHOT2(4)
      REAL*8  R,X(4),Y(4),PI0(4)
      INTEGER JEZELI(3),K
      DATA JEZELI /0,0,0/
      SAVE JEZELI

! random 3 vector on the sphere, masless
        R=SQRT(PI0(4)**2-PI0(3)**2-PI0(2)**2-PI0(1)**2)/2D0
        CALL SPHERD(R,X)
        X(4)=R
        Y(4)=R
        
        Y(1)=-X(1)
        Y(2)=-X(2)
        Y(3)=-X(3)
! boost to lab and to real*4
        CALL bostdq(-1,PI0,X,X)
        CALL bostdq(-1,PI0,Y,Y)
        DO L=1,4
         PHOT1(L)=X(L)
         PHOT2(L)=Y(L)
        ENDDO
C to hepevt
        CALL FILHEP(0,1,22,K,K,0,0,PHOT1,0.0,.TRUE.)
        CALL FILHEP(0,1,22,K,K,0,0,PHOT2,0.0,.TRUE.)

C
      END
      SUBROUTINE TAUETA(PETA,K)
C subroutine to decay etas's from taus. 
C this routine belongs to tauola universal interface, but uses 
C routines from tauola utilities. Just flat phase space, but 4 channels.
C it is called at the beginning of SUBR. TAUPI0(JAK)
C and as far as hepevt search it is basically the same as TAUPI0.  25.08.2005    
C this is the hepevt class in old style. No d_h_ class pre-name
 
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

C position of taus, must be defined by host program:
c
      REAL  RRR(1),BRSUM(3), RR(2)
      REAL  PHOT1(4),PHOT2(4),PHOT3(4)
      REAL*8    X(4),    Y(4),    Z(4)
      REAL                                YM1,YM2,YM3
      REAL*8  R,RU,PETA(4),XM1,XM2,XM3,XLAM
      REAL*8 a,b,c
      INTEGER K
      XLAM(a,b,c)=SQRT(ABS((a-b-c)**2-4.0*b*c))
C position of decaying particle:
C        DO L=1,4
C          PETA(L)= phep(L,K)  ! eta 4 momentum
C        ENDDO
C       eta cumulated branching ratios:
        BRSUM(1)=0.389  ! gamma gamma
        BRSUM(2)=BRSUM(1)+0.319  ! 3 pi0
        BRSUM(3)=BRSUM(2)+0.237  ! pi+ pi- pi0 rest is thus pi+pi-gamma
        CALL RANMAR(RRR,1) 
        
        IF (RRR(1).LT.BRSUM(1)) THEN ! gamma gamma channel exactly like pi0
! random 3 vector on the sphere, masless   
         R=SQRT(PETA(4)**2-PETA(3)**2-PETA(2)**2-PETA(1)**2)/2D0
         CALL SPHERD(R,X) 
         X(4)=R
         Y(4)=R
        
         Y(1)=-X(1)
         Y(2)=-X(2)
         Y(3)=-X(3)
! boost to lab and to real*4
         CALL bostdq(-1,PETA,X,X)  
         CALL bostdq(-1,PETA,Y,Y)
         DO L=1,4
          PHOT1(L)=X(L)
          PHOT2(L)=Y(L)
         ENDDO
C to hepevt
         CALL FILHEP(0,1,22,K,K,0,0,PHOT1,0.0,.TRUE.)
         CALL FILHEP(0,1,22,K,K,0,0,PHOT2,0.0,.TRUE.)
        ELSE ! 3 body channels
         IF(RRR(1).LT.BRSUM(2)) THEN  ! 3 pi0
          ID1= 111
          ID2= 111
          ID3= 111
          XM1=AMPIZ ! masses
          XM2=AMPIZ
          XM3=AMPIZ
         ELSEIF(RRR(1).LT.BRSUM(3)) THEN ! pi+ pi- pi0
          ID1= 211
          ID2=-211
          ID3= 111
          XM1=AMPI ! masses
          XM2=AMPI
          XM3=AMPIZ
         ELSE                            ! pi+ pi- gamma 
          ID1= 211
          ID2=-211
          ID3=  22
          XM1=AMPI ! masses
          XM2=AMPI
          XM3=0.0
         ENDIF
 7       CONTINUE  ! we generate mass of the first pair:
          CALL RANMAR(RR,2)
          R=SQRT(PETA(4)**2-PETA(3)**2-PETA(2)**2-PETA(1)**2)
          AMIN=XM1+XM2
          AMAX=R-XM3
          AM2=SQRT(AMIN**2+RR(1)*(AMAX**2-AMIN**2))
C         weight for flat phase space
          WT=XLAM(1D0*R**2,1D0*AM2**2,1D0*XM3**2)
     &      *XLAM(1D0*AM2**2,1D0*XM1**2,1D0*XM2**2)
     &           /R**2                    /AM2**2
         IF (RR(2).GT.WT) GOTO 7

         RU=XLAM(1D0*AM2**2,1D0*XM1**2,1D0*XM2**2)/AM2/2  ! momenta of the
                                              ! first two products
                                              ! in the rest frame of that pair
         CALL SPHERD(RU,X)
         X(4)=SQRT(RU**2+XM1**2)
         Y(4)=SQRT(RU**2+XM2**2)
        
         Y(1)=-X(1)
         Y(2)=-X(2)
         Y(3)=-X(3)
C generate momentum of that pair in rest frame of eta:
         RU=XLAM(1D0*R**2,1D0*AM2**2,1D0*XM3**2)/R/2
         CALL SPHERD(RU,Z)
         Z(4)=SQRT(RU**2+AM2**2)
C and boost first two decay products to rest frame of eta.
         CALL bostdq(-1,Z,X,X)
         CALL bostdq(-1,Z,Y,Y)
C redefine Z(4) to 4-momentum of the last decay product: 
         Z(1)=-Z(1)
         Z(2)=-Z(2)
         Z(3)=-Z(3)
         Z(4)=SQRT(RU**2+XM3**2)
C boost all to lab and move to real*4; also masses
         CALL bostdq(-1,PETA,X,X)
         CALL bostdq(-1,PETA,Y,Y)
         CALL bostdq(-1,PETA,Z,Z)
         DO L=1,4
          PHOT1(L)=X(L)
          PHOT2(L)=Y(L)
          PHOT3(L)=Z(L)
         ENDDO
         YM1=XM1
         YM2=XM2
         YM3=XM3
C to hepevt
         CALL FILHEP(0,1,ID1,K,K,0,0,PHOT1,YM1,.TRUE.)
         CALL FILHEP(0,1,ID2,K,K,0,0,PHOT2,YM2,.TRUE.)
         CALL FILHEP(0,1,ID3,K,K,0,0,PHOT3,YM3,.TRUE.)
        ENDIF


C
      END
      SUBROUTINE TAUK0S(PETA,K)
C subroutine to decay K0S's from taus. 
C this routine belongs to tauola universal interface, but uses 
C routines from tauola utilities. Just flat phase space, but 4 channels.
C it is called at the beginning of SUBR. TAUPI0(JAK)
C and as far as hepevt search it is basically the same as TAUPI0.  25.08.2005   

      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU 
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST

C position of taus, must be defined by host program:
      COMMON /TAUPOS/ NP1,NP2
c
      REAL  RRR(1),BRSUM(3)
      REAL  PHOT1(4),PHOT2(4)
      REAL*8    X(4),    Y(4)
      REAL                                YM1,YM2
      REAL*8  R,PETA(4),XM1,XM2,XLAM
      REAL*8 a,b,c
      INTEGER K
      XLAM(a,b,c)=SQRT(ABS((a-b-c)**2-4.0*b*c))
C position of decaying particle:

      
!        DO L=1,4
!          PETA(L)= phep(L,K)  ! K0S 4 momentum  (this is cloned from eta decay)
!        ENDDO
C       K0S cumulated branching ratios:
        BRSUM(1)=0.313  ! 2 PI0
        BRSUM(2)=1.0 ! BRSUM(1)+0.319  ! Pi+ PI-
        BRSUM(3)=BRSUM(2)+0.237  ! pi+ pi- pi0 rest is thus pi+pi-gamma
        CALL RANMAR(RRR,1) 

         IF(RRR(1).LT.BRSUM(1)) THEN  ! 2 pi0
          ID1= 111
          ID2= 111
          XM1=AMPIZ ! masses
          XM2=AMPIZ
         ELSEIF(RRR(1).LT.BRSUM(2)) THEN ! pi+ pi- 
          ID1= 211
          ID2=-211
          XM1=AMPI ! masses
          XM2=AMPI
         ELSE                            ! gamma gamma unused !!!
          ID1= 22
          ID2= 22
          XM1= 0.0 ! masses
          XM2= 0.0
         ENDIF
        
! random 3 vector on the sphere, of equal mass !!  
         R=SQRT(PETA(4)**2-PETA(3)**2-PETA(2)**2-PETA(1)**2)/2D0
         R4=R
         R=SQRT(ABS(R**2-XM1**2))
         CALL SPHERD(R,X) 
         X(4)=R4
         Y(4)=R4
        
         Y(1)=-X(1)
         Y(2)=-X(2)
         Y(3)=-X(3)
! boost to lab and to real*4
         CALL bostdq(-1,PETA,X,X)  
         CALL bostdq(-1,PETA,Y,Y)
         DO L=1,4
          PHOT1(L)=X(L)
          PHOT2(L)=Y(L)
         ENDDO

         YM1=XM1
         YM2=XM2
C to hepevt
         CALL FILHEP(0,1,ID1,K,K,0,0,PHOT1,YM1,.TRUE.)
         CALL FILHEP(0,1,ID2,K,K,0,0,PHOT2,YM2,.TRUE.)

C
      END

      subroutine bostdq(idir,vv,pp,q)
*     *******************************
c Boost along arbitrary vector v (see eg. J.D. Jacson, Classical 
c Electrodynamics).
c Four-vector pp is boosted from an actual frame to the rest frame 
c of the four-vector v (for idir=1) or back (for idir=-1). 
c q is a resulting four-vector.
c Note: v must be time-like, pp may be arbitrary.
c
c Written by: Wieslaw Placzek            date: 22.07.1994
c Last update: 3/29/95                     by: M.S.
c 
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (nout=6)
      DOUBLE PRECISION v(4),p(4),q(4),pp(4),vv(4)  
      save
!
      do 1 i=1,4
      v(i)=vv(i)
 1    p(i)=pp(i)
      amv=(v(4)**2-v(1)**2-v(2)**2-v(3)**2)
      if (amv.le.0d0) then
        write(6,*) 'bosstv: warning amv**2=',amv
      endif
      amv=sqrt(abs(amv))
      if (idir.eq.-1) then
        q(4)=( p(1)*v(1)+p(2)*v(2)+p(3)*v(3)+p(4)*v(4))/amv
        wsp =(q(4)+p(4))/(v(4)+amv)
      elseif (idir.eq.1) then
        q(4)=(-p(1)*v(1)-p(2)*v(2)-p(3)*v(3)+p(4)*v(4))/amv
        wsp =-(q(4)+p(4))/(v(4)+amv)
      else
        write(nout,*)' >>> boostv: wrong value of idir = ',idir
      endif
      q(1)=p(1)+wsp*v(1)
      q(2)=p(2)+wsp*v(2)
      q(3)=p(3)+wsp*v(3)
      end
        




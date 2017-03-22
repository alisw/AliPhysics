* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYSSPA_HIJING(IPU1,IPU2)  
    
C...Generates spacelike parton showers. 
      IMPLICIT DOUBLE PRECISION(D)  
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint3_hijing.inc"
      DIMENSION KFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVS(2),ROBO(5),   
     &XFS(2,-6:6),XFA(-6:6),XFB(-6:6),XFN(-6:6),WTAP(-6:6),WTSF(-6:6),  
     &THE2(2),ALAM(2),DQ2(3),DPC(3),DPD(4),DPB(4)   
    
C...Calculate maximum virtuality and check that evolution possible. 
      IPUS1=IPU1    
      IPUS2=IPU2    
      ISUB=MINT(1)  
      Q2E=VINT(52)  
      IF(ISET(ISUB).EQ.1) THEN  
        Q2E=Q2E/PARP(67)    
      ELSEIF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN   
        Q2E=PMAS(23,1)**2   
        IF(ISUB.EQ.8.OR.ISUB.EQ.76.OR.ISUB.EQ.77) Q2E=PMAS(24,1)**2 
      ENDIF 
      TMAX=LOG(PARP(67)*PARP(63)*Q2E/PARP(61)**2)   
      IF(PARP(67)*Q2E.LT.MAX(PARP(62)**2,2.*PARP(61)**2).OR.    
     &TMAX.LT.0.2) RETURN   
    
C...Common constants and initial values. Save normal Lambda value.  
      XE0=2.*PARP(65)/VINT(1)   
      ALAMS=PARU(111)   
      PARU(111)=PARP(61)    
      NS=N  
  100 N=NS  
      DO 110 JT=1,2 
      KFLS(JT)=MINT(14+JT)  
      KFLS(JT+2)=KFLS(JT)   
      XS(JT)=VINT(40+JT)    
      ZS(JT)=1. 
      Q2S(JT)=PARP(67)*Q2E  
      TEVS(JT)=TMAX 
      ALAM(JT)=PARP(61) 
      THE2(JT)=100. 
      DO 110 KFL=-6,6   
  110 XFS(JT,KFL)=XSFX(JT,KFL)  
      DSH=VINT(44)  
      IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) DSH=VINT(26)*VINT(2)   
    
C...Pick up leg with highest virtuality.    
  120 N=N+1 
      JT=1  
      IF(N.GT.NS+1.AND.Q2S(2).GT.Q2S(1)) JT=2   
      KFLB=KFLS(JT) 
      XB=XS(JT) 
      DO 130 KFL=-6,6   
  130 XFB(KFL)=XFS(JT,KFL)  
      DSHR=2D0*SQRT(DSH)    
      DSHZ=DSH/DBLE(ZS(JT)) 
      XE=MAX(XE0,XB*(1./(1.-PARP(66))-1.))  
      IF(XB+XE.GE.0.999) THEN   
        Q2B=0.  
        GOTO 220    
      ENDIF 
    
C...Maximum Q2 without or with Q2 ordering. Effective Lambda and n_f.   
      IF(MSTP(62).LE.1) THEN    
        Q2B=0.5*(1./ZS(JT)+1.)*Q2S(JT)+0.5*(1./ZS(JT)-1.)*(Q2S(3-JT)-   
     &  SNGL(DSH)+SQRT((SNGL(DSH)+Q2S(1)+Q2S(2))**2+8.*Q2S(1)*Q2S(2)*   
     &  ZS(JT)/(1.-ZS(JT))))    
        TEVB=LOG(PARP(63)*Q2B/ALAM(JT)**2)  
      ELSE  
        Q2B=Q2S(JT) 
        TEVB=TEVS(JT)   
      ENDIF 
      ALSDUM=ULALPS_HIJING(PARP(63)*Q2B)   
      TEVB=TEVB+2.*LOG(ALAM(JT)/PARU(117))  
      TEVBSV=TEVB   
      ALAM(JT)=PARU(117)    
      B0=(33.-2.*MSTU(118))/6.  
    
C...Calculate Altarelli-Parisi and structure function weights.  
      DO 140 KFL=-6,6   
      WTAP(KFL)=0.  
  140 WTSF(KFL)=0.  
      IF(KFLB.EQ.21) THEN   
        WTAPQ=16.*(1.-SQRT(XB+XE))/(3.*SQRT(XB))    
        DO 150 KFL=-MSTP(54),MSTP(54)   
        IF(KFL.EQ.0) WTAP(KFL)=6.*LOG((1.-XB)/XE)   
  150   IF(KFL.NE.0) WTAP(KFL)=WTAPQ    
      ELSE  
        WTAP(0)=0.5*XB*(1./(XB+XE)-1.)  
        WTAP(KFLB)=8.*LOG((1.-XB)*(XB+XE)/XE)/3.    
      ENDIF 
  160 WTSUM=0.  
      IF(KFLB.NE.21) XFBO=XFB(KFLB) 
      IF(KFLB.EQ.21) XFBO=XFB(0)
C***************************************************************
C**********ERROR HAS OCCURED HERE
      IF(XFBO.EQ.0.0) THEN
                WRITE(MSTU(11),1000)
                WRITE(MSTU(11),1001) KFLB,XFB(KFLB)
		XFBO=0.00001
      ENDIF
C****************************************************************    
      DO 170 KFL=-MSTP(54),MSTP(54) 
      WTSF(KFL)=XFB(KFL)/XFBO   
  170 WTSUM=WTSUM+WTAP(KFL)*WTSF(KFL)   
      WTSUM=MAX(0.0001,WTSUM)   
    
C...Choose new t: fix alpha_s, alpha_s(Q2), alpha_s(k_T2).  
  180 IF(MSTP(64).LE.0) THEN    
        TEVB=TEVB+LOG(RLU_HIJING(0))*PARU(2)/(PARU(111)*WTSUM) 
      ELSEIF(MSTP(64).EQ.1) THEN    
        TEVB=TEVB*EXP(MAX(-100.,LOG(RLU_HIJING(0))*B0/WTSUM))  
      ELSE  
        TEVB=TEVB*EXP(MAX(-100.,LOG(RLU_HIJING(0))*B0/(5.*WTSUM))) 
      ENDIF 
  190 Q2REF=ALAM(JT)**2*EXP(TEVB)   
      Q2B=Q2REF/PARP(63)    
    
C...Evolution ended or select flavour for branching parton. 
      IF(Q2B.LT.PARP(62)**2) THEN   
        Q2B=0.  
      ELSE  
        WTRAN=RLU_HIJING(0)*WTSUM  
        KFLA=-MSTP(54)-1    
  200   KFLA=KFLA+1 
        WTRAN=WTRAN-WTAP(KFLA)*WTSF(KFLA)   
        IF(KFLA.LT.MSTP(54).AND.WTRAN.GT.0.) GOTO 200   
        IF(KFLA.EQ.0) KFLA=21   
    
C...Choose z value and corrective weight.   
        IF(KFLB.EQ.21.AND.KFLA.EQ.21) THEN  
          Z=1./(1.+((1.-XB)/XB)*(XE/(1.-XB))**RLU_HIJING(0))   
          WTZ=(1.-Z*(1.-Z))**2  
        ELSEIF(KFLB.EQ.21) THEN 
          Z=XB/(1.-RLU_HIJING(0)*(1.-SQRT(XB+XE)))**2  
          WTZ=0.5*(1.+(1.-Z)**2)*SQRT(Z)    
        ELSEIF(KFLA.EQ.21) THEN 
          Z=XB*(1.+RLU_HIJING(0)*(1./(XB+XE)-1.))  
          WTZ=1.-2.*Z*(1.-Z)    
        ELSE    
          Z=1.-(1.-XB)*(XE/((XB+XE)*(1.-XB)))**RLU_HIJING(0)   
          WTZ=0.5*(1.+Z**2) 
        ENDIF   
    
C...Option with resummation of soft gluon emission as effective z shift.    
        IF(MSTP(65).GE.1) THEN  
          RSOFT=6.  
          IF(KFLB.NE.21) RSOFT=8./3.    
          Z=Z*(TEVB/TEVS(JT))**(RSOFT*XE/((XB+XE)*B0))  
          IF(Z.LE.XB) GOTO 180  
        ENDIF   
    
C...Option with alpha_s(k_T2)Q2): demand k_T2 > cutoff, reweight.   
        IF(MSTP(64).GE.2) THEN  
          IF((1.-Z)*Q2B.LT.PARP(62)**2) GOTO 180    
          ALPRAT=TEVB/(TEVB+LOG(1.-Z))  
          IF(ALPRAT.LT.5.*RLU_HIJING(0)) GOTO 180  
          IF(ALPRAT.GT.5.) WTZ=WTZ*ALPRAT/5.    
        ENDIF   
    
C...Option with angular ordering requirement.   
        IF(MSTP(62).GE.3) THEN  
          THE2T=(4.*Z**2*Q2B)/(VINT(2)*(1.-Z)*XB**2)    
          IF(THE2T.GT.THE2(JT)) GOTO 180    
        ENDIF   
    
C...Weighting with new structure functions. 
        CALL PYSTFU_HIJING(MINT(10+JT),XB,Q2REF,XFN,JT)   
        IF(KFLB.NE.21) XFBN=XFN(KFLB)   
        IF(KFLB.EQ.21) XFBN=XFN(0)  
        IF(XFBN.LT.1E-20) THEN  
          IF(KFLA.EQ.KFLB) THEN 
            TEVB=TEVBSV 
            WTAP(KFLB)=0.   
            GOTO 160    
          ELSEIF(TEVBSV-TEVB.GT.0.2) THEN   
            TEVB=0.5*(TEVBSV+TEVB)  
            GOTO 190    
          ELSE  
            XFBN=1E-10  
          ENDIF 
        ENDIF   
        DO 210 KFL=-MSTP(54),MSTP(54)   
  210   XFB(KFL)=XFN(KFL)   
        XA=XB/Z 
        CALL PYSTFU_HIJING(MINT(10+JT),XA,Q2REF,XFA,JT)   
        IF(KFLA.NE.21) XFAN=XFA(KFLA)   
        IF(KFLA.EQ.21) XFAN=XFA(0)  
        IF(XFAN.LT.1E-20) GOTO 160  
        IF(KFLA.NE.21) WTSFA=WTSF(KFLA) 
        IF(KFLA.EQ.21) WTSFA=WTSF(0)    
        IF(WTZ*XFAN/XFBN.LT.RLU_HIJING(0)*WTSFA) GOTO 160  
      ENDIF 
    
C...Define two hard scatterers in their CM-frame.   
  220 IF(N.EQ.NS+2) THEN    
        DQ2(JT)=Q2B 
        DPLCM=SQRT((DSH+DQ2(1)+DQ2(2))**2-4D0*DQ2(1)*DQ2(2))/DSHR   
        DO 240 JR=1,2   
        I=NS+JR 
        IF(JR.EQ.1) IPO=IPUS1   
        IF(JR.EQ.2) IPO=IPUS2   
        DO 230 J=1,5    
        K(I,J)=0    
        P(I,J)=0.   
  230   V(I,J)=0.   
        K(I,1)=14   
        K(I,2)=KFLS(JR+2)   
        K(I,4)=IPO  
        K(I,5)=IPO  
        P(I,3)=DPLCM*(-1)**(JR+1)   
        P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR 
        P(I,5)=-SQRT(SNGL(DQ2(JR))) 
        K(IPO,1)=14 
        K(IPO,3)=I  
        K(IPO,4)=MOD(K(IPO,4),MSTU(5))+MSTU(5)*I    
  240   K(IPO,5)=MOD(K(IPO,5),MSTU(5))+MSTU(5)*I    
    
C...Find maximum allowed mass of timelike parton.   
      ELSEIF(N.GT.NS+2) THEN    
        JR=3-JT 
        DQ2(3)=Q2B  
        DPC(1)=P(IS(1),4)   
        DPC(2)=P(IS(2),4)   
        DPC(3)=0.5*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))    
        DPD(1)=DSH+DQ2(JR)+DQ2(JT)  
        DPD(2)=DSHZ+DQ2(JR)+DQ2(3)  
        DPD(3)=SQRT(DPD(1)**2-4D0*DQ2(JR)*DQ2(JT))  
        DPD(4)=SQRT(DPD(2)**2-4D0*DQ2(JR)*DQ2(3))   
        IKIN=0  
        IF(Q2S(JR).GE.(0.5*PARP(62))**2.AND.DPD(1)-DPD(3).GE.   
     &  1D-10*DPD(1)) IKIN=1    
        IF(IKIN.EQ.0) DMSMA=(DQ2(JT)/DBLE(ZS(JT))-DQ2(3))*(DSH/ 
     &  (DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))    
        IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/(2.*  
     &  DQ2(JR))-DQ2(JT)-DQ2(3) 
    
C...Generate timelike parton shower (if required).  
        IT=N    
        DO 250 J=1,5    
        K(IT,J)=0   
        P(IT,J)=0.  
  250   V(IT,J)=0.  
        K(IT,1)=3   
        K(IT,2)=21  
        IF(KFLB.EQ.21.AND.KFLS(JT+2).NE.21) K(IT,2)=-KFLS(JT+2) 
        IF(KFLB.NE.21.AND.KFLS(JT+2).EQ.21) K(IT,2)=KFLB    
        P(IT,5)=ULMASS_HIJING(K(IT,2)) 
        IF(SNGL(DMSMA).LE.P(IT,5)**2) GOTO 100  
        IF(MSTP(63).GE.1) THEN  
          P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR    
          P(IT,3)=SQRT(P(IT,4)**2-P(IT,5)**2)   
          IF(MSTP(63).EQ.1) THEN    
            Q2TIM=DMSMA 
          ELSEIF(MSTP(63).EQ.2) THEN    
            Q2TIM=MIN(SNGL(DMSMA),PARP(71)*Q2S(JT)) 
          ELSE  
C'''Here remains to introduce angular ordering in first branching.  
            Q2TIM=DMSMA 
          ENDIF 
          CALL LUSHOW_HIJING(IT,0,SQRT(Q2TIM)) 
          IF(N.GE.IT+1) P(IT,5)=P(IT+1,5)   
        ENDIF   
    
C...Reconstruct kinematics of branching: timelike parton shower.    
        DMS=P(IT,5)**2  
        IF(IKIN.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/(DSH+DQ2(JT))  
        IF(IKIN.EQ.1) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*DPD(2)+0.5*DPD(3)*   
     &  DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/(4.*DSH*DPC(3)**2) 
        IF(DPT2.LT.0.) GOTO 100 
        DPB(1)=(0.5*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/  
     &  DSHR)/DPC(3)-DPC(3) 
        P(IT,1)=SQRT(SNGL(DPT2))    
        P(IT,3)=DPB(1)*(-1)**(JT+1) 
        P(IT,4)=(DSHZ-DSH-DMS)/DSHR 
        IF(N.GE.IT+1) THEN  
          DPB(1)=SQRT(DPB(1)**2+DPT2)   
          DPB(2)=SQRT(DPB(1)**2+DMS)    
          DPB(3)=P(IT+1,3)  
          DPB(4)=SQRT(DPB(3)**2+DMS)    
          DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)* 
     &    DPB(1))   
          CALL LUDBRB_HIJING(IT+1,N,0.,0.,0D0,0D0,DBEZ)    
          THE=ULANGL_HIJING(P(IT,3),P(IT,1))   
          CALL LUDBRB_HIJING(IT+1,N,THE,0.,0D0,0D0,0D0)    
        ENDIF   
    
C...Reconstruct kinematics of branching: spacelike parton.  
        DO 260 J=1,5    
        K(N+1,J)=0  
        P(N+1,J)=0. 
  260   V(N+1,J)=0. 
        K(N+1,1)=14 
        K(N+1,2)=KFLB   
        P(N+1,1)=P(IT,1)    
        P(N+1,3)=P(IT,3)+P(IS(JT),3)    
        P(N+1,4)=P(IT,4)+P(IS(JT),4)    
        P(N+1,5)=-SQRT(SNGL(DQ2(3)))    
    
C...Define colour flow of branching.    
        K(IS(JT),3)=N+1 
        K(IT,3)=N+1 
        ID1=IT  
        IF((K(N+1,2).GT.0.AND.K(N+1,2).NE.21.AND.K(ID1,2).GT.0.AND. 
     &  K(ID1,2).NE.21).OR.(K(N+1,2).LT.0.AND.K(ID1,2).EQ.21).OR.   
     &  (K(N+1,2).EQ.21.AND.K(ID1,2).EQ.21.AND.RLU_HIJING(0).GT.0.5).OR.   
     &  (K(N+1,2).EQ.21.AND.K(ID1,2).LT.0)) ID1=IS(JT)  
        ID2=IT+IS(JT)-ID1   
        K(N+1,4)=K(N+1,4)+ID1   
        K(N+1,5)=K(N+1,5)+ID2   
        K(ID1,4)=K(ID1,4)+MSTU(5)*(N+1) 
        K(ID1,5)=K(ID1,5)+MSTU(5)*ID2   
        K(ID2,4)=K(ID2,4)+MSTU(5)*ID1   
        K(ID2,5)=K(ID2,5)+MSTU(5)*(N+1) 
        N=N+1   
    
C...Boost to new CM-frame.  
        CALL LUDBRB_HIJING(NS+1,N,0.,0.,-DBLE((P(N,1)+P(IS(JR),1))/(P(N
     $       ,4)+P(IS(JR),4))),0D0,-DBLE((P(N,3)+P(IS(JR),3))/(P(N,4)
     $       +P(IS(JR),4))))  
        IR=N+(JT-1)*(IS(1)-N)   
        CALL LUDBRB_HIJING(NS+1,N,-ULANGL_HIJING(P(IR,3),P(IR,1)),PARU(2
     $       )*RLU_HIJING(0),0D0,0D0,0D0)    
      ENDIF 
    
C...Save quantities, loop back. 
      IS(JT)=N  
      Q2S(JT)=Q2B   
      DQ2(JT)=Q2B   
      IF(MSTP(62).GE.3) THE2(JT)=THE2T  
      DSH=DSHZ  
      IF(Q2B.GE.(0.5*PARP(62))**2) THEN 
        KFLS(JT+2)=KFLS(JT) 
        KFLS(JT)=KFLA   
        XS(JT)=XA   
        ZS(JT)=Z    
        DO 270 KFL=-6,6 
  270   XFS(JT,KFL)=XFA(KFL)    
        TEVS(JT)=TEVB   
      ELSE  
        IF(JT.EQ.1) IPU1=N  
        IF(JT.EQ.2) IPU2=N  
      ENDIF 
      IF(N.GT.MSTU(4)-MSTU(32)-10) THEN 
         CALL LUERRM_HIJING(11
     $        ,'(PYSSPA_HIJING:) no more memory left in LUJETS_HIJING')   
        IF(MSTU(21).GE.1) N=NS  
        IF(MSTU(21).GE.1) RETURN    
      ENDIF 
      IF(MAX(Q2S(1),Q2S(2)).GE.(0.5*PARP(62))**2.OR.N.LE.NS+1) GOTO 120 
    
C...Boost hard scattering partons to frame of shower initiators.    
      DO 280 J=1,3  
  280 ROBO(J+2)=(P(NS+1,J)+P(NS+2,J))/(P(NS+1,4)+P(NS+2,4)) 
      DO 290 J=1,5  
  290 P(N+2,J)=P(NS+1,J)    
      ROBOT=ROBO(3)**2+ROBO(4)**2+ROBO(5)**2    
      IF(ROBOT.GE.0.999999) THEN    
        ROBOT=1.00001*SQRT(ROBOT)   
        ROBO(3)=ROBO(3)/ROBOT   
        ROBO(4)=ROBO(4)/ROBOT   
        ROBO(5)=ROBO(5)/ROBOT   
      ENDIF 
      CALL LUDBRB_HIJING(N+2,N+2,0.,0.,-DBLE(ROBO(3)),-DBLE(ROBO(4)),  
     &-DBLE(ROBO(5)))   
      ROBO(2)=ULANGL_HIJING(P(N+2,1),P(N+2,2)) 
      ROBO(1)=ULANGL_HIJING(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))    
      CALL LUDBRB_HIJING(MINT(83)+5,NS,ROBO(1),ROBO(2),DBLE(ROBO(3)),  
     &DBLE(ROBO(4)),DBLE(ROBO(5)))  
    
C...Store user information. Reset Lambda value. 
      K(IPU1,3)=MINT(83)+3  
      K(IPU2,3)=MINT(83)+4  
      DO 300 JT=1,2 
      MINT(12+JT)=KFLS(JT)  
  300 VINT(140+JT)=XS(JT)   
      PARU(111)=ALAMS   
 1000 FORMAT(5X,'structure function has a zero point here')
 1001 FORMAT(5X,'xf(x,i=',I5,')=',F10.5)

      RETURN    
      END   

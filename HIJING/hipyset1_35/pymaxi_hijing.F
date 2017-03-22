* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYMAXI_HIJING 
    
C...Finds optimal set of coefficients for kinematical variable selection    
C...and the maximum of the part of the differential cross-section used  
C...in the event weighting. 
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "pysubs_hijing.inc"
#include "pypars_hijing.inc"
#include "pyint1_hijing.inc"
#include "pyint2_hijing.inc"
#include "pyint3_hijing.inc"
#include "pyint4_hijing.inc"
#include "pyint5_hijing.inc"
#include "pyint6_hijing.inc"
      CHARACTER CVAR(4)*4   
      DIMENSION NPTS(4),MVARPT(200,4),VINTPT(200,30),SIGSPT(200),   
     &NAREL(6),WTREL(6),WTMAT(6,6),COEFU(6),IACCMX(4),SIGSMX(4),    
     &SIGSSM(3) 
      DATA CVAR/'tau ','tau''','y*  ','cth '/   
    
C...Select subprocess to study: skip cases not applicable.  
      VINT(143)=1.  
      VINT(144)=1.  
      XSEC(0,1)=0.  
      DO 350 ISUB=1,200 
      IF(ISUB.GE.91.AND.ISUB.LE.95) THEN    
        XSEC(ISUB,1)=VINT(ISUB+11)  
        IF(MSUB(ISUB).NE.1) GOTO 350    
        GOTO 340    
      ELSEIF(ISUB.EQ.96) THEN   
        IF(MINT(43).NE.4) GOTO 350  
        IF(MSUB(95).NE.1.AND.MSTP(81).LE.0.AND.MSTP(131).LE.0) GOTO 350 
      ELSEIF(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.ISUB.EQ.28.OR.   
     &ISUB.EQ.53.OR.ISUB.EQ.68) THEN    
        IF(MSUB(ISUB).NE.1.OR.MSUB(95).EQ.1) GOTO 350   
      ELSE  
        IF(MSUB(ISUB).NE.1) GOTO 350    
      ENDIF 
      MINT(1)=ISUB  
      ISTSB=ISET(ISUB)  
      IF(ISUB.EQ.96) ISTSB=2    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1000) ISUB  
    
C...Find resonances (explicit or implicit in cross-section).    
      MINT(72)=0    
      KFR1=0    
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3) THEN 
        KFR1=KFPR(ISUB,1)   
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN    
        KFR1=25 
      ENDIF 
      IF(KFR1.NE.0) THEN    
        TAUR1=PMAS(KFR1,1)**2/VINT(2)   
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2) 
        MINT(72)=1  
        MINT(73)=KFR1   
        VINT(73)=TAUR1  
        VINT(74)=GAMR1  
      ENDIF 
      IF(ISUB.EQ.141) THEN  
        KFR2=23 
        TAUR2=PMAS(KFR2,1)**2/VINT(2)   
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2) 
        MINT(72)=2  
        MINT(74)=KFR2   
        VINT(75)=TAUR2  
        VINT(76)=GAMR2  
      ENDIF 
    
C...Find product masses and minimum pT of process.  
      SQM3=0.   
      SQM4=0.   
      MINT(71)=0    
      VINT(71)=CKIN(3)  
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN 
        IF(KFPR(ISUB,1).NE.0) SQM3=PMAS(KFPR(ISUB,1),1)**2  
        IF(KFPR(ISUB,2).NE.0) SQM4=PMAS(KFPR(ISUB,2),1)**2  
        IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) MINT(71)=1 
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5)) 
        IF(ISUB.EQ.96.AND.MSTP(82).LE.1) VINT(71)=PARP(81)  
        IF(ISUB.EQ.96.AND.MSTP(82).GE.2) VINT(71)=0.08*PARP(82) 
      ENDIF 
      VINT(63)=SQM3 
      VINT(64)=SQM4 
    
C...Number of points for each variable: tau, tau', y*, cos(theta-hat).  
      NPTS(1)=2+2*MINT(72)  
      IF(MINT(43).EQ.1.AND.(ISTSB.EQ.1.OR.ISTSB.EQ.2)) NPTS(1)=1    
      NPTS(2)=1 
      IF(MINT(43).GE.2.AND.(ISTSB.EQ.3.OR.ISTSB.EQ.4)) NPTS(2)=2    
      NPTS(3)=1 
      IF(MINT(43).EQ.4) NPTS(3)=3   
      NPTS(4)=1 
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) NPTS(4)=5    
      NTRY=NPTS(1)*NPTS(2)*NPTS(3)*NPTS(4)  
    
C...Reset coefficients of cross-section weighting.  
      DO 100 J=1,20 
  100 COEF(ISUB,J)=0.   
      COEF(ISUB,1)=1.   
      COEF(ISUB,7)=0.5  
      COEF(ISUB,8)=0.5  
      COEF(ISUB,10)=1.  
      COEF(ISUB,15)=1.  
      MCTH=0    
      MTAUP=0   
      CTH=0.    
      TAUP=0.   
      SIGSAM=0. 
    
C...Find limits and select tau, y*, cos(theta-hat) and tau' values, 
C...in grid of phase space points.  
      CALL PYKLIM_HIJING(1)    
      NACC=0    
      DO 120 ITRY=1,NTRY    
      IF(MOD(ITRY-1,NPTS(2)*NPTS(3)*NPTS(4)).EQ.0) THEN 
        MTAU=1+(ITRY-1)/(NPTS(2)*NPTS(3)*NPTS(4))   
        CALL PYKMAP_HIJING(1,MTAU,0.5) 
        IF(ISTSB.EQ.3.OR.ISTSB.EQ.4) CALL PYKLIM_HIJING(4) 
      ENDIF 
      IF((ISTSB.EQ.3.OR.ISTSB.EQ.4).AND.MOD(ITRY-1,NPTS(3)*NPTS(4)).    
     &EQ.0) THEN    
        MTAUP=1+MOD((ITRY-1)/(NPTS(3)*NPTS(4)),NPTS(2)) 
        CALL PYKMAP_HIJING(4,MTAUP,0.5)    
      ENDIF 
      IF(MOD(ITRY-1,NPTS(3)*NPTS(4)).EQ.0) CALL PYKLIM_HIJING(2)   
      IF(MOD(ITRY-1,NPTS(4)).EQ.0) THEN 
        MYST=1+MOD((ITRY-1)/NPTS(4),NPTS(3))    
        CALL PYKMAP_HIJING(2,MYST,0.5) 
        CALL PYKLIM_HIJING(3)  
      ENDIF 
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN 
        MCTH=1+MOD(ITRY-1,NPTS(4))  
        CALL PYKMAP_HIJING(3,MCTH,0.5) 
      ENDIF 
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2) 
    
C...Calculate and store cross-section.  
      MINT(51)=0    
      CALL PYKLIM_HIJING(0)    
      IF(MINT(51).EQ.1) GOTO 120    
      NACC=NACC+1   
      MVARPT(NACC,1)=MTAU   
      MVARPT(NACC,2)=MTAUP  
      MVARPT(NACC,3)=MYST   
      MVARPT(NACC,4)=MCTH   
      DO 110 J=1,30 
  110 VINTPT(NACC,J)=VINT(10+J) 
      CALL PYSIGH_HIJING(NCHN,SIGS)    
      SIGSPT(NACC)=SIGS 
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1100) MTAU,MTAUP,MYST,MCTH, 
     &VINT(21),VINT(22),VINT(23),VINT(26),SIGS  
  120 CONTINUE  
      IF(SIGSAM.EQ.0.) THEN 
        WRITE(MSTU(11),1200) ISUB   
        STOP    
      ENDIF 
    
C...Calculate integrals in tau and y* over maximal phase space limits.  
      TAUMIN=VINT(11)   
      TAUMAX=VINT(31)   
      ATAU1=LOG(TAUMAX/TAUMIN)  
      ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN) 
      IF(NPTS(1).GE.3) THEN 
        ATAU3=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1    
        ATAU4=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/  
     &  GAMR1   
      ENDIF 
      IF(NPTS(1).GE.5) THEN 
        ATAU5=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2    
        ATAU6=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/  
     &  GAMR2   
      ENDIF 
      YSTMIN=0.5*LOG(TAUMIN)    
      YSTMAX=-YSTMIN    
      AYST0=YSTMAX-YSTMIN   
      AYST1=0.5*(YSTMAX-YSTMIN)**2  
      AYST3=2.*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))    
    
C...Reset. Sum up cross-sections in points calculated.  
      DO 230 IVAR=1,4   
      IF(NPTS(IVAR).EQ.1) GOTO 230  
      IF(ISUB.EQ.96.AND.IVAR.EQ.4) GOTO 230 
      NBIN=NPTS(IVAR)   
      DO 130 J1=1,NBIN  
      NAREL(J1)=0   
      WTREL(J1)=0.  
      COEFU(J1)=0.  
      DO 130 J2=1,NBIN  
  130 WTMAT(J1,J2)=0.   
      DO 140 IACC=1,NACC    
      IBIN=MVARPT(IACC,IVAR)    
      NAREL(IBIN)=NAREL(IBIN)+1 
      WTREL(IBIN)=WTREL(IBIN)+SIGSPT(IACC)  
    
C...Sum up tau cross-section pieces in points used. 
      IF(IVAR.EQ.1) THEN    
        TAU=VINTPT(IACC,11) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAU1/ATAU2)/TAU   
        IF(NBIN.GE.3) THEN  
          WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAU1/ATAU3)/(TAU+TAUR1) 
          WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ATAU1/ATAU4)*TAU/    
     &    ((TAU-TAUR1)**2+GAMR1**2) 
        ENDIF   
        IF(NBIN.GE.5) THEN  
          WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ATAU1/ATAU5)/(TAU+TAUR2) 
          WTMAT(IBIN,6)=WTMAT(IBIN,6)+(ATAU1/ATAU6)*TAU/    
     &    ((TAU-TAUR2)**2+GAMR2**2) 
        ENDIF   
    
C...Sum up tau' cross-section pieces in points used.    
      ELSEIF(IVAR.EQ.2) THEN    
        TAU=VINTPT(IACC,11) 
        TAUP=VINTPT(IACC,16)    
        TAUPMN=VINTPT(IACC,6)   
        TAUPMX=VINTPT(IACC,26)  
        ATAUP1=LOG(TAUPMX/TAUPMN)   
        ATAUP2=((1.-TAU/TAUPMX)**4-(1.-TAU/TAUPMN)**4)/(4.*TAU) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAUP1/ATAUP2)*(1.-TAU/TAUP)**3/   
     &  TAUP    
    
C...Sum up y* and cos(theta-hat) cross-section pieces in points used.   
      ELSEIF(IVAR.EQ.3) THEN    
        YST=VINTPT(IACC,12) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+(AYST0/AYST1)*(YST-YSTMIN)  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(AYST0/AYST1)*(YSTMAX-YST)  
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(AYST0/AYST3)/COSH(YST) 
      ELSE  
        RM34=2.*SQM3*SQM4/(VINTPT(IACC,11)*VINT(2))**2  
        RSQM=1.+RM34    
        CTHMAX=SQRT(1.-4.*VINT(71)**2/(TAUMAX*VINT(2))) 
        CTHMIN=-CTHMAX  
        IF(CTHMAX.GT.0.9999) RM34=MAX(RM34,2.*VINT(71)**2/  
     &  (TAUMAX*VINT(2)))   
        ACTH1=CTHMAX-CTHMIN 
        ACTH2=LOG(MAX(RM34,RSQM-CTHMIN)/MAX(RM34,RSQM-CTHMAX))  
        ACTH3=LOG(MAX(RM34,RSQM+CTHMAX)/MAX(RM34,RSQM+CTHMIN))  
        ACTH4=1./MAX(RM34,RSQM-CTHMAX)-1./MAX(RM34,RSQM-CTHMIN) 
        ACTH5=1./MAX(RM34,RSQM+CTHMIN)-1./MAX(RM34,RSQM+CTHMAX) 
        CTH=VINTPT(IACC,13) 
        WTMAT(IBIN,1)=WTMAT(IBIN,1)+1.  
        WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ACTH1/ACTH2)/MAX(RM34,RSQM-CTH)    
        WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ACTH1/ACTH3)/MAX(RM34,RSQM+CTH)    
        WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ACTH1/ACTH4)/MAX(RM34,RSQM-CTH)**2 
        WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ACTH1/ACTH5)/MAX(RM34,RSQM+CTH)**2 
      ENDIF 
  140 CONTINUE  
    
C...Check that equation system solvable; else trivial way out.  
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1300) CVAR(IVAR)    
      MSOLV=1   
      DO 150 IBIN=1,NBIN    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1400) (WTMAT(IBIN,IRED),    
     &IRED=1,NBIN),WTREL(IBIN)  
  150 IF(NAREL(IBIN).EQ.0) MSOLV=0  
      IF(MSOLV.EQ.0) THEN   
        DO 160 IBIN=1,NBIN  
  160   COEFU(IBIN)=1.  
    
C...Solve to find relative importance of cross-section pieces.  
      ELSE  
        DO 170 IRED=1,NBIN-1    
        DO 170 IBIN=IRED+1,NBIN 
        RQT=WTMAT(IBIN,IRED)/WTMAT(IRED,IRED)   
        WTREL(IBIN)=WTREL(IBIN)-RQT*WTREL(IRED) 
        DO 170 ICOE=IRED,NBIN   
  170   WTMAT(IBIN,ICOE)=WTMAT(IBIN,ICOE)-RQT*WTMAT(IRED,ICOE)  
        DO 190 IRED=NBIN,1,-1   
        DO 180 ICOE=IRED+1,NBIN 
  180   WTREL(IRED)=WTREL(IRED)-WTMAT(IRED,ICOE)*COEFU(ICOE)    
  190   COEFU(IRED)=WTREL(IRED)/WTMAT(IRED,IRED)    
      ENDIF 
    
C...Normalize coefficients, with piece shared democratically.   
      COEFSU=0. 
      DO 200 IBIN=1,NBIN    
      COEFU(IBIN)=MAX(0.,COEFU(IBIN))   
  200 COEFSU=COEFSU+COEFU(IBIN) 
      IF(IVAR.EQ.1) IOFF=0  
      IF(IVAR.EQ.2) IOFF=14 
      IF(IVAR.EQ.3) IOFF=6  
      IF(IVAR.EQ.4) IOFF=9  
      IF(COEFSU.GT.0.) THEN 
        DO 210 IBIN=1,NBIN  
  210   COEF(ISUB,IOFF+IBIN)=PARP(121)/NBIN+(1.-PARP(121))*COEFU(IBIN)/ 
     &  COEFSU  
      ELSE  
        DO 220 IBIN=1,NBIN  
  220   COEF(ISUB,IOFF+IBIN)=1./NBIN    
      ENDIF 
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1500) CVAR(IVAR),   
     &(COEF(ISUB,IOFF+IBIN),IBIN=1,NBIN)    
  230 CONTINUE  
    
C...Find two most promising maxima among points previously determined.  
      DO 240 J=1,4  
      IACCMX(J)=0   
  240 SIGSMX(J)=0.  
      NMAX=0    
      DO 290 IACC=1,NACC    
      DO 250 J=1,30 
  250 VINT(10+J)=VINTPT(IACC,J) 
      CALL PYSIGH_HIJING(NCHN,SIGS)    
      IEQ=0 
      DO 260 IMV=1,NMAX 
  260 IF(ABS(SIGS-SIGSMX(IMV)).LT.1E-4*(SIGS+SIGSMX(IMV))) IEQ=IMV  
      IF(IEQ.EQ.0) THEN 
        DO 270 IMV=NMAX,1,-1    
        IIN=IMV+1   
        IF(SIGS.LE.SIGSMX(IMV)) GOTO 280    
        IACCMX(IMV+1)=IACCMX(IMV)   
  270   SIGSMX(IMV+1)=SIGSMX(IMV)   
        IIN=1   
  280   IACCMX(IIN)=IACC    
        SIGSMX(IIN)=SIGS    
        IF(NMAX.LE.1) NMAX=NMAX+1   
      ENDIF 
  290 CONTINUE  
    
C...Read out starting position for search.  
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1600)   
      SIGSAM=SIGSMX(1)  
      DO 330 IMAX=1,NMAX    
      IACC=IACCMX(IMAX) 
      MTAU=MVARPT(IACC,1)   
      MTAUP=MVARPT(IACC,2)  
      MYST=MVARPT(IACC,3)   
      MCTH=MVARPT(IACC,4)   
      VTAU=0.5  
      VYST=0.5  
      VCTH=0.5  
      VTAUP=0.5 
    
C...Starting point and step size in parameter space.    
      DO 320 IRPT=1,2   
      DO 310 IVAR=1,4   
      IF(NPTS(IVAR).EQ.1) GOTO 310  
      IF(IVAR.EQ.1) VVAR=VTAU   
      IF(IVAR.EQ.2) VVAR=VTAUP  
      IF(IVAR.EQ.3) VVAR=VYST   
      IF(IVAR.EQ.4) VVAR=VCTH   
      IF(IVAR.EQ.1) MVAR=MTAU   
      IF(IVAR.EQ.2) MVAR=MTAUP  
      IF(IVAR.EQ.3) MVAR=MYST   
      IF(IVAR.EQ.4) MVAR=MCTH   
      IF(IRPT.EQ.1) VDEL=0.1    
      IF(IRPT.EQ.2) VDEL=MAX(0.01,MIN(0.05,VVAR-0.02,0.98-VVAR))    
      IF(IRPT.EQ.1) VMAR=0.02   
      IF(IRPT.EQ.2) VMAR=0.002  
      IMOV0=1   
      IF(IRPT.EQ.1.AND.IVAR.EQ.1) IMOV0=0   
      DO 300 IMOV=IMOV0,8   
    
C...Define new point in parameter space.    
      IF(IMOV.EQ.0) THEN    
        INEW=2  
        VNEW=VVAR   
      ELSEIF(IMOV.EQ.1) THEN    
        INEW=3  
        VNEW=VVAR+VDEL  
      ELSEIF(IMOV.EQ.2) THEN    
        INEW=1  
        VNEW=VVAR-VDEL  
      ELSEIF(SIGSSM(3).GE.MAX(SIGSSM(1),SIGSSM(2)).AND. 
     &VVAR+2.*VDEL.LT.1.-VMAR) THEN 
        VVAR=VVAR+VDEL  
        SIGSSM(1)=SIGSSM(2) 
        SIGSSM(2)=SIGSSM(3) 
        INEW=3  
        VNEW=VVAR+VDEL  
      ELSEIF(SIGSSM(1).GE.MAX(SIGSSM(2),SIGSSM(3)).AND. 
     &VVAR-2.*VDEL.GT.VMAR) THEN    
        VVAR=VVAR-VDEL  
        SIGSSM(3)=SIGSSM(2) 
        SIGSSM(2)=SIGSSM(1) 
        INEW=1  
        VNEW=VVAR-VDEL  
      ELSEIF(SIGSSM(3).GE.SIGSSM(1)) THEN   
        VDEL=0.5*VDEL   
        VVAR=VVAR+VDEL  
        SIGSSM(1)=SIGSSM(2) 
        INEW=2  
        VNEW=VVAR   
      ELSE  
        VDEL=0.5*VDEL   
        VVAR=VVAR-VDEL  
        SIGSSM(3)=SIGSSM(2) 
        INEW=2  
        VNEW=VVAR   
      ENDIF 
    
C...Convert to relevant variables and find derived new limits.  
      IF(IVAR.EQ.1) THEN    
        VTAU=VNEW   
        CALL PYKMAP_HIJING(1,MTAU,VTAU)    
        IF(ISTSB.EQ.3.OR.ISTSB.EQ.4) CALL PYKLIM_HIJING(4) 
      ENDIF 
      IF(IVAR.LE.2.AND.(ISTSB.EQ.3.OR.ISTSB.EQ.4)) THEN 
        IF(IVAR.EQ.2) VTAUP=VNEW    
        CALL PYKMAP_HIJING(4,MTAUP,VTAUP)  
      ENDIF 
      IF(IVAR.LE.2) CALL PYKLIM_HIJING(2)  
      IF(IVAR.LE.3) THEN    
        IF(IVAR.EQ.3) VYST=VNEW 
        CALL PYKMAP_HIJING(2,MYST,VYST)    
        CALL PYKLIM_HIJING(3)  
      ENDIF 
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN 
        IF(IVAR.EQ.4) VCTH=VNEW 
        CALL PYKMAP_HIJING(3,MCTH,VCTH)    
      ENDIF 
      IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2) 
    
C...Evaluate cross-section. Save new maximum. Final maximum.    
      CALL PYSIGH_HIJING(NCHN,SIGS)    
      SIGSSM(INEW)=SIGS 
      IF(SIGS.GT.SIGSAM) SIGSAM=SIGS    
      IF(MSTP(122).GE.2) WRITE(MSTU(11),1700) IMAX,IVAR,MVAR,IMOV,  
     &VNEW,VINT(21),VINT(22),VINT(23),VINT(26),SIGS 
  300 CONTINUE  
  310 CONTINUE  
  320 CONTINUE  
      IF(IMAX.EQ.1) SIGS11=SIGSAM   
  330 CONTINUE  
      XSEC(ISUB,1)=1.05*SIGSAM  
  340 IF(ISUB.NE.96) XSEC(0,1)=XSEC(0,1)+XSEC(ISUB,1)   
  350 CONTINUE  
    
C...Print summary table.    
      IF(MSTP(122).GE.1) THEN   
        WRITE(MSTU(11),1800)    
        WRITE(MSTU(11),1900)    
        DO 360 ISUB=1,200   
        IF(MSUB(ISUB).NE.1.AND.ISUB.NE.96) GOTO 360 
        IF(ISUB.EQ.96.AND.MINT(43).NE.4) GOTO 360   
        IF(ISUB.EQ.96.AND.MSUB(95).NE.1.AND.MSTP(81).LE.0) GOTO 360 
        IF(MSUB(95).EQ.1.AND.(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.    
     &  ISUB.EQ.28.OR.ISUB.EQ.53.OR.ISUB.EQ.68)) GOTO 360   
        WRITE(MSTU(11),2000) ISUB,PROC(ISUB),XSEC(ISUB,1)   
  360   CONTINUE    
        WRITE(MSTU(11),2100)    
      ENDIF 
    
C...Format statements for maximization results. 
 1000 FORMAT(/1X,'Coefficient optimization and maximum search for ',    
     &'subprocess no',I4/1X,'Coefficient modes     tau',10X,'y*',9X,    
     &'cth',9X,'tau''',7X,'sigma')  
 1100 FORMAT(1X,4I4,F12.8,F12.6,F12.7,F12.8,1P,E12.4)   
 1200 FORMAT(1X,'Error: requested subprocess ',I3,' has vanishing ',    
     &'cross-section.'/1X,'Execution stopped!')
 1300 FORMAT(1X,'Coefficients of equation system to be solved for ',A4) 
 1400 FORMAT(1X,1P,7E11.3)  
 1500 FORMAT(1X,'Result for ',A4,':',6F9.4) 
 1600 FORMAT(1X,'Maximum search for given coefficients'/2X,'MAX VAR ',  
     &'MOD MOV   VNEW',7X,'tau',7X,'y*',8X,'cth',7X,'tau''',7X,'sigma') 
 1700 FORMAT(1X,4I4,F8.4,F11.7,F9.3,F11.6,F11.7,1P,E12.4)   
 1800 FORMAT(/1X,8('*'),1X,'PYMAXI_HIJING: summary of differential ',  
     &'cross-section maximum search',1X,8('*')) 
 1900 FORMAT(/11X,58('=')/11X,'I',38X,'I',17X,'I'/11X,'I  ISUB  ',  
     &'Subprocess name',15X,'I  Maximum value  I'/11X,'I',38X,'I',  
     &17X,'I'/11X,58('=')/11X,'I',38X,'I',17X,'I')  
 2000 FORMAT(11X,'I',2X,I3,3X,A28,2X,'I',2X,1P,E12.4,3X,'I')    
 2100 FORMAT(11X,'I',38X,'I',17X,'I'/11X,58('='))   
    
      RETURN    
      END   

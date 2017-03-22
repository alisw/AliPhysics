* $Id$
    
C*********************************************************************  
    
      SUBROUTINE PYSTFE_HIJING(KF,X,Q2,XPQ)    
    
C...This is a dummy routine, where the user can introduce an interface  
C...to his own external structure function parametrization. 
C...Arguments in:   
C...KF : 2212 for p, 211 for pi+; isospin conjugation for n and charge  
C...    conjugation for pbar, nbar or pi- is performed by PYSTFU.   
C...X : x value.    
C...Q2 : Q^2 value. 
C...Arguments out:  
C...XPQ(-6:6) : x * f(x,Q2), with index according to KF code,   
C...    except that gluon is placed in 0. Thus XPQ(0) = xg, 
C...    XPQ(1) = xd, XPQ(-1) = xdbar, XPQ(2) = xu, XPQ(-2) = xubar, 
C...    XPQ(3) = xs, XPQ(-3) = xsbar, XPQ(4) = xc, XPQ(-4) = xcbar, 
C...    XPQ(5) = xb, XPQ(-5) = xbbar, XPQ(6) = xt, XPQ(-6) = xtbar. 
C...    
C...One such interface, to the Diemos, Ferroni, Longo, Martinelli   
C...proton structure functions, already comes with the package. What    
C...the user needs here is external files with the three routines   
C...FXG160, FXG260 and FXG360 of the authors above, plus the    
C...interpolation routine FINT, which is part of the CERN library   
C...KERNLIB package. To avoid problems with unresolved external 
C...references, the external calls are commented in the current 
C...version. To enable this option, remove the C* at the beginning  
C...of the relevant lines.  
C...    
C...Alternatively, the routine can be used as an interface to the   
C...structure function evolution program of Tung. This can be achieved  
C...by removing C* at the beginning of some of the lines below. 
#include "ludat1_hijing.inc"
#include "ludat2_hijing.inc"
#include "pypars_hijing.inc"
      DIMENSION XPQ(-6:6),XFDFLM(9) 
      CHARACTER CHDFLM(9)*5,HEADER*40   
      DATA CHDFLM/'UPVAL','DOVAL','GLUON','QBAR ','UBAR ','SBAR ',  
     &'CBAR ','BBAR ','TBAR '/  
      DATA HEADER/'Tung evolution package has been invoked'/    
      DATA INIT/0/  
    
C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli. 
C...Allowed variable range 10 GeV2 < Q2 < 1E8 GeV2, 5E-5 < x < .95. 
      IF(MSTP(51).GE.11.AND.MSTP(51).LE.13.AND.MSTP(52).LE.1) THEN  
        XDFLM=MAX(0.51E-4,X)    
        Q2DFLM=MAX(10.,MIN(1E8,Q2)) 
        IF(MSTP(52).EQ.0) Q2DFLM=10.    
        DO 100 J=1,9    
        IF(MSTP(52).EQ.1.AND.J.EQ.9) THEN   
          Q2DFLM=Q2DFLM*(40./PMAS(6,1))**2  
          Q2DFLM=MAX(10.,MIN(1E8,Q2))   
        ENDIF   
        XFDFLM(J)=0.    
C...Remove C* on following three lines to enable the DFLM options.  
C*      IF(MSTP(51).EQ.11) CALL FXG160(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))    
C*      IF(MSTP(51).EQ.12) CALL FXG260(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))    
C*      IF(MSTP(51).EQ.13) CALL FXG360(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))    
  100   CONTINUE    
        IF(X.LT.0.51E-4.AND.ABS(PARP(51)-1.).GT.0.01) THEN  
          CXS=(0.51E-4/X)**(PARP(51)-1.)    
          DO 110 J=1,7  
  110     XFDFLM(J)=XFDFLM(J)*CXS   
        ENDIF   
        XPQ(0)=XFDFLM(3)    
        XPQ(1)=XFDFLM(2)+XFDFLM(5)  
        XPQ(2)=XFDFLM(1)+XFDFLM(5)  
        XPQ(3)=XFDFLM(6)    
        XPQ(4)=XFDFLM(7)    
        XPQ(5)=XFDFLM(8)    
        XPQ(6)=XFDFLM(9)    
        XPQ(-1)=XFDFLM(5)   
        XPQ(-2)=XFDFLM(5)   
        XPQ(-3)=XFDFLM(6)   
        XPQ(-4)=XFDFLM(7)   
        XPQ(-5)=XFDFLM(8)   
        XPQ(-6)=XFDFLM(9)   
    
C...Proton structure function evolution from Wu-Ki Tung: parton 
C...distribution functions incorporating heavy quark mass effects.  
C...Allowed variable range: PARP(52) < Q < PARP(53); PARP(54) < x < 1.  
      ELSE  
        IF(INIT.EQ.0) THEN  
          I1=0  
          IF(MSTP(52).EQ.4) I1=1    
          IHDRN=1   
          NU=MSTP(53)   
          I2=MSTP(51)   
          IF(MSTP(51).GE.11) I2=MSTP(51)-3  
          I3=0  
          IF(MSTP(52).EQ.3) I3=1    
    
C...Convert to Lambda in CWZ scheme (approximately linear relation).    
          ALAM=0.75*PARP(1) 
          TPMS=PMAS(6,1)    
          QINI=PARP(52) 
          QMAX=PARP(53) 
          XMIN=PARP(54) 
    
C...Initialize evolution (perform calculation or read results from  
C...file).  
C...Remove C* on following two lines to enable Tung initialization. 
C*        CALL PDFSET(I1,IHDRN,ALAM,TPMS,QINI,QMAX,XMIN,NU,HEADER,  
C*   &    I2,I3,IRET,IRR)   
          INIT=1    
        ENDIF   
    
C...Put into output array.  
        Q=SQRT(Q2)  
        DO 200 I=-6,6   
        FIXQ=0. 
C...Remove C* on following line to enable structure function call.  
C*      FIXQ=MAX(0.,PDF(10,1,I,X,Q,IR)) 
  200   XPQ(I)=X*FIXQ   
    
C...Change order of u and d quarks from Tung to PYTHIA convention.  
        XPS=XPQ(1)  
        XPQ(1)=XPQ(2)   
        XPQ(2)=XPS  
        XPS=XPQ(-1) 
        XPQ(-1)=XPQ(-2) 
        XPQ(-2)=XPS 
      ENDIF 
    
      RETURN    
      END   

#include "isajet/pilot.h"
      LOGICAL FUNCTION DECJET(IP,NADD,IDABS,PREST,WDECAY,BETA,GAMMA)
C
C          Auxiliary routine for DECAY. Evolve and hadronize partons.
C          Check conservation laws. Return TRUE if OK, FALSE otherwise.
C
C          IP = particle to be decayed.
C          NADD = number of products (NPTCL+1, ..., NPTCL+NADD).
C          IDABS = absolute values of decay IDENT's.
C          PREST = 4-momenta in rest frame.
C          WDECAY = logical flag for real W decay.
C          BETA,GAMMA = boost parameters.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/wcon.inc"
#include "isajet/partcl.inc"
#include "isajet/dkytab.inc"
#include "isajet/jetset.inc"
#include "isajet/jwork.inc"
#include "isajet/const.inc"
C   
      REAL PGEN(5,5),RND(5),U(3),BETA(3),IDQK(3),ROT(3,3),PSAVE(3)  
     1,REDUCE(5),WPROP,Z,TRY,RANF,AMASS,TWOME,CHARGE    
      REAL PSUM(5),POLD(4),PNEW(4),SUM,WTMAX,SUM1,SUM2  
      REAL PREST(4,6),PWREST(5),BETAW(3),DOT,PCM    
      REAL AMEE,REE,WTEE,SWAP,RNEW,WT,QCM,PHI,S12,S12MAX,GAMMAW,BP  
      REAL PJET,CTHQK,STHQK,CPHIQK,SPHIQK,SUMQ,A,B,C,GAMMA  
      REAL CHARGW   
      LOGICAL WDECAY    
      INTEGER IDLV1,IFL1,IFL2,IFL3,JSPIN,INDEX,IPOINT,ID1,I1,I2,I3  
      INTEGER NADD,NSTART,NEW,NADD1,J,IP,I,IDABS(5),NEXT    
      INTEGER JJ1,II,K1,K,NJSAVE,NJSAV1,NJSAV2,NJ1,NPRTN,KK,NHDRN1  
      INTEGER IFAIL,JSAVE,JETIP,JET,NJADD,NPTLV1,IDANTI,NPJET(5)    
      INTEGER NHDRN,NPJET3,NPTCLW,NPBEG(5)  
C   
C          Copy decay products into /JETSET/ and do QCD evolution.  
C   
      IF(NJSET+NADD.GT.MXJSET) GO TO 9998   
      NJSAVE=NJSET  
      NSTART=NPTCL-NADD+1   
      NPTCL=NSTART-1    
      DO 100 I=1,NADD   
        NJSET=NJSET+1   
        DO 110 K=1,4    
110     PJSET(K,NJSET)=PREST(K,I)   
        PJSET(5,NJSET)=PPTCL(5,NPTCL+I) 
        JORIG(NJSET)=JPACK*I    
        JTYPE(NJSET)=IDENT(NPTCL+I) 
        JDCAY(NJSET)=0  
        JMATCH(NJSET)=JPACK*(NJSAVE+1)+NJSAVE+NADD  
100   CONTINUE  
C   
C          For heavy quarks match 1+2 and 3+(1+2). Boost 1+2 to rest.   
C   
      IF(WDECAY) THEN   
        JMATCH(NJSAVE+1)=NJSAVE+2   
        JMATCH(NJSAVE+2)=NJSAVE+1   
        NJSET=NJSET+1   
        DO 120 K=1,4    
          PWREST(K)=PJSET(K,NJSAVE+1)+PJSET(K,NJSAVE+2) 
          PJSET(K,NJSET)=PWREST(K)  
120     CONTINUE    
        PWREST(5)=SQRT(PWREST(4)**2-PWREST(1)**2-PWREST(2)**2   
     $  -PWREST(3)**2)  
        PJSET(5,NJSET)=PWREST(5)    
        JMATCH(NJSAVE+3)=NJSAVE+4   
        JMATCH(NJSAVE+4)=NJSAVE+3   
        JORIG(NJSAVE+4)=-1  
        IDLV1=JTYPE(NJSAVE+1)   
        CHARGW=CHARGE(IDLV1)    
        IDLV1=JTYPE(NJSAVE+2)   
        CHARGW=CHARGW+CHARGE(IDLV1) 
        JTYPE(NJSAVE+4)=80*SIGN(1.,CHARGW)  
        JDCAY(NJSAVE+4)=0   
C          Boost W vectors to rest. 
        DO 130 K=1,3    
130     BETAW(K)=PWREST(K)/PWREST(4)    
        GAMMAW=PWREST(4)/PWREST(5)  
        NJSAV1=NJSAVE+1 
        NJSAV2=NJSAVE+2 
        DO 140 J=NJSAV1,NJSAV2  
          BP=BETAW(1)*PJSET(1,J)+BETAW(2)*PJSET(2,J)+BETAW(3)*PJSET(3,J)    
          DO 141 K=1,3  
141       PJSET(K,J)=PJSET(K,J)-GAMMAW*BETAW(K)*(PJSET(4,J) 
     $    -BP*GAMMAW/(GAMMAW+1.))   
          PJSET(4,J)=GAMMAW*(PJSET(4,J)-BP) 
140     CONTINUE    
      ENDIF 
C   
C          Do evolution and save new W momentum. Start from parent  
C          mass or NADD*energy. 
      NJSAV1=NJSAVE+1   
      DO 150 J=NJSAV1,NJSET 
        IF(IABS(JTYPE(J)).LT.10.OR.MOD(JTYPE(J),100).EQ.0) THEN 
          JDCAY(J)=-1   
          PJSET(5,J)=AMIN1(PPTCL(5,IP),NADD*PJSET(4,J)) 
        ENDIF   
150   CONTINUE  
C   
      CALL QCDJET(NJSAVE+1) 
C   
      IF(WDECAY) THEN   
        PWREST(4)=PJSET(4,NJSAVE+4) 
        GAMMAW=PWREST(4)/PWREST(5)  
        DO 200 K=1,3    
          PWREST(K)=PJSET(K,NJSAVE+4)   
          BETAW(K)=PWREST(K)/PWREST(4)  
200     CONTINUE    
      ENDIF 
C   
C          Put final partons in particle table - temporary IORIG.   
C          Also include virtual or real W momentum for quark decays.    
C   
      NJ1=NJSAVE+1  
      IF(WDECAY) THEN   
C          Real or virtual W.   
        NPTCL=NPTCL+1   
        NPTCLW=NPTCL    
        DO 210 K=1,5    
210     PPTCL(K,NPTCL)=PJSET(K,NJSAVE+4)    
        IORIG(NPTCL)=IP 
        IDENT(NPTCL)=JTYPE(NJSAVE+4)    
        IDCAY(NPTCL)=0  
C          Jet 3    
        NPBEG(3)=NPTCL+1    
        DO 220 J=NJ1,NJSET  
          IF(JDCAY(J).NE.0) GO TO 220   
          IF(JORIG(J)/JPACK.NE.3) GO TO 220 
          NPTCL=NPTCL+1 
          DO 221 K=1,5  
221       PPTCL(K,NPTCL)=PJSET(K,J) 
          IORIG(NPTCL)=3*IPACK+IP   
          IDENT(NPTCL)=JTYPE(J) 
          IDCAY(NPTCL)=0    
220     CONTINUE    
C          Jets 1 and 2 
        NPJET3=NPTCL    
        DO 230 JET=1,2  
          NPBEG(JET)=NPTCL+1    
          DO 240 J=NJ1,NJSET    
            IF(JDCAY(J).NE.0) GO TO 240 
            IF(JORIG(J)/JPACK.NE.JET) GO TO 240 
            NPTCL=NPTCL+1   
            BP=BETAW(1)*PJSET(1,J)+BETAW(2)*PJSET(2,J)  
     $      +BETAW(3)*PJSET(3,J)    
            DO 241 K=1,3    
241         PPTCL(K,NPTCL)=PJSET(K,J)+GAMMAW*BETAW(K)*(PJSET(4,J)   
     $      +BP*GAMMAW/(GAMMAW+1.)) 
            PPTCL(4,NPTCL)=GAMMAW*(PJSET(4,J)+BP)   
            PPTCL(5,NPTCL)=PJSET(5,J)   
            IORIG(NPTCL)=IPACK*(JORIG(J)/JPACK)+NPTCLW  
            IDENT(NPTCL)=JTYPE(J)   
            IDCAY(NPTCL)=0  
240       CONTINUE  
230     CONTINUE    
C          Quark decays to W plus jet 3; then W decays. 
        IDCAY(IP)=IPACK*NPTCLW+NPJET3   
        IDCAY(NPTCLW)=IPACK*(NPJET3+1)+NPTCL    
      ELSE  
C          Not quark decay, so just copy partons.   
        DO 250 JET=1,NADD   
          NPBEG(JET)=NPTCL+1    
          DO 260 J=NJ1,NJSET    
            IF(JDCAY(J).NE.0) GO TO 260 
            IF(JORIG(J)/JPACK.NE.JET) GO TO 260 
            NPTCL=NPTCL+1   
            DO 261 K=1,5    
261         PPTCL(K,NPTCL)=PJSET(K,J)   
            IORIG(NPTCL)=IPACK*(JORIG(J)/JPACK)+IP  
            IDENT(NPTCL)=JTYPE(J)   
            IDCAY(NPTCL)=0  
260       CONTINUE  
250     CONTINUE    
        IDCAY(IP)=NSTART*IPACK+NPTCL    
      ENDIF 
      NHDRN=NPTCL   
C   
C          Hadronize quarks and rotate to proper angles.    
C   
      DO 300 JET=1,NADD 
        NPRTN=NPBEG(JET)-1  
        DO 310 I=NJ1,NJSET  
          IF(JDCAY(I).NE.0) GO TO 310   
          IF(JORIG(I)/JPACK.NE.JET) GO TO 310   
          NPRTN=NPRTN+1 
          IF(IABS(JTYPE(I)).GE.10.AND.MOD(JTYPE(I),100).NE.0)   
     $    GO TO 330 
C   
C          Fragment parton: 
          NEXT=NPTCL+1  
          PJET=SQRT(PJSET(1,I)**2+PJSET(2,I)**2+PJSET(3,I)**2)  
          CTHQK=PJSET(3,I)/PJET 
          STHQK=1.-CTHQK**2
          IF(STHQK.LT.1) THEN
            STHQK=SQRT(STHQK)
            CPHIQK=PJSET(1,I)/(PJET*STHQK)
            SPHIQK=PJSET(2,I)/(PJET*STHQK)
          ELSE
            STHQK=0
            CPHIQK=1
            SPHIQK=0
          ENDIF
          CALL JETGEN(I)    
          IF(NEXT.GT.NPTCL) GO TO 310   
          ROT(1,1)=CPHIQK*CTHQK 
          ROT(2,1)=SPHIQK*CTHQK 
          ROT(3,1)=-STHQK   
          ROT(1,2)=-SPHIQK  
          ROT(2,2)=CPHIQK   
          ROT(3,2)=0.   
          ROT(1,3)=CPHIQK*STHQK 
          ROT(2,3)=SPHIQK*STHQK 
          ROT(3,3)=CTHQK    
C   
          DO 320 II=NEXT,NPTCL  
            DO 321 K=1,3    
              PSAVE(K)=PPTCL(K,II)  
              PPTCL(K,II)=0.    
321         CONTINUE    
            DO 322 K=1,3    
            DO 322 KK=1,3   
322         PPTCL(K,II)=PPTCL(K,II)+ROT(K,KK)*PSAVE(KK) 
            IORIG(II)=IPACK*JET+NPRTN   
            IDCAY(II)=0 
320       CONTINUE  
          IDCAY(NPRTN)=NEXT*IPACK+NPTCL 
          GO TO 310 
C   
C          or add lepton:   
330       NPTCL=NPTCL+1 
          DO 331 K=1,5  
331       PPTCL(K,NPTCL)=PJSET(K,I) 
          IORIG(NPTCL)=IPACK*JET+NPRTN  
          IDENT(NPTCL)=JTYPE(I) 
          IDCAY(NPTCL)=0    
          IDCAY(NPRTN)=NPTCL*IPACK+NPTCL    
310     CONTINUE    
        NPJET(JET)=NPTCL    
300   CONTINUE  
C   
C          Reset NJSET so decay jets do not appear in /JETSET/  
      NJADD=NJSET   
      NJSET=NJSAVE  
C   
C          Check for at least two particles 
      IF(NPTCL.LT.NHDRN+2) THEN 
        NPTCL=NSTART-1  
        DECJET=.FALSE.  
        RETURN  
      ENDIF 
C   
C          Conserve charge  
C   
      SUMQ=0.   
      NHDRN1=NHDRN+1    
      DO 400 I=NHDRN1,NPTCL 
        IDLV1=IDENT(I)  
        SUMQ=SUMQ+CHARGE(IDLV1) 
400   CONTINUE  
      IDLV1=IDENT(IP)   
      SUMQ=SUMQ-CHARGE(IDLV1)   
C   
      IF(SUMQ.EQ.0.) GO TO 500  
C   
C          Charge wrong--fix it by swapping UP and DN quarks.   
      DO 410 I=NHDRN1,NPTCL 
        ID1=IDENT(I)    
        IF(IABS(ID1).GT.1000) GO TO 410 
        I1=MOD(IABS(ID1)/100,10)    
        I2=MOD(IABS(ID1)/10,10) 
        I3=MOD(IABS(ID1),10)    
        IF(I1.EQ.1.AND.I2.GT.2.AND.SUMQ*ID1.GT.0.) THEN 
          IDENT(I)=ISIGN(200+10*I2+I3,ID1)  
        ELSEIF(I1.EQ.2.AND.I2.GT.2.AND.SUMQ*ID1.LT.0.) THEN 
          IDENT(I)=ISIGN(100+10*I2+I3,ID1)  
        ELSEIF(I1.EQ.1.AND.I2.EQ.2.AND.SUMQ*ID1.GT.0.) THEN 
          IDENT(I)=110+I3   
        ELSEIF(I1.EQ.1.AND.I2.EQ.1) THEN    
          IDENT(I)=(120+I3)*(-SIGN(1.,SUMQ))    
        ELSE    
          GO TO 410 
        ENDIF   
        SUMQ=SIGN(ABS(SUMQ)-1.,SUMQ)    
        IDLV1=IDENT(I)  
        PPTCL(5,I)=AMASS(IDLV1) 
        PPTCL(4,I)=SQRT(PPTCL(1,I)**2+PPTCL(2,I)**2+PPTCL(3,I)**2   
     $  +PPTCL(5,I)**2) 
C          Sum cannot vanish for fractionally charged initial particle. 
        IF(ABS(SUMQ).LT.1.) GO TO 500   
410   CONTINUE  
C          Failed to conserve charge.   
      NPTCL=NSTART-1    
      DECJET=.FALSE.    
      RETURN    
C   
C          Rescale momenta for correct mass 
C   
500   CONTINUE  
      IF(WDECAY) THEN   
C          Quark decay. First rescale jet3 + W  
        DO 510 K=1,5    
510     PPTCL(K,NPTCL+1)=PPTCL(K,NPTCLW)    
        NPTLV1=NPTCL+1  
        DO 520 K=1,3    
520     PSUM(K)=0.  
        PSUM(4)=PPTCL(5,IP) 
        PSUM(5)=PSUM(4) 
        CALL RESCAL(NPJET(2)+1,NPTLV1,PSUM,IFAIL) 
        IF(IFAIL.NE.0) THEN 
          NPTCL=NSTART-1    
          DECJET=.FALSE.    
          RETURN    
        ENDIF   
        DO 530 K=1,3    
530     BETAW(K)=PPTCL(K,NPTCL+1)/PPTCL(4,NPTCL+1)  
        GAMMAW=PPTCL(4,NPTCL+1)/PPTCL(5,NPTCL+1)    
C          Then rescale W   
        PSUM(4)=PPTCL(5,NPTCLW) 
        PSUM(5)=PSUM(4) 
        CALL RESCAL(NHDRN1,NPJET(2),PSUM,IFAIL) 
        IF(IFAIL.NE.0) THEN 
          NPTCL=NSTART-1    
          DECJET=.FALSE.    
          RETURN    
        ENDIF   
      ELSE  
C          General decay with no W. 
        DO 550 K=1,3    
550     PSUM(K)=0.  
        PSUM(4)=PPTCL(5,IP) 
        PSUM(5)=PSUM(4) 
        NPTLV1=NPTCL    
        CALL RESCAL(NHDRN1,NPTLV1,PSUM,IFAIL)   
        IF(IFAIL.NE.0) THEN 
          NPTCL=NSTART-1    
          DECJET=.FALSE.    
          RETURN    
        ENDIF   
      ENDIF 
C   
C          Boost back to lab frame. Reset IORIG.    
C   
      IF(WDECAY) THEN   
        DO 600 I=NHDRN1,NPTCL  
          JET=IORIG(I)/IPACK    
          IF(JET.NE.1.AND.JET.NE.2) GO TO 600   
          BP=BETAW(1)*PPTCL(1,I)+BETAW(2)*PPTCL(2,I)+BETAW(3)*PPTCL(3,I)    
          DO 610 J=1,3  
610       PPTCL(J,I)=PPTCL(J,I)+GAMMAW*BETAW(J)*(PPTCL(4,I) 
     $    +BP*GAMMAW/(GAMMAW+1.))   
          PPTCL(4,I)=GAMMAW*(PPTCL(4,I)+BP) 
600     CONTINUE    
      ENDIF 
C   
      DO 620 I=NSTART,NPTCL 
        IORIG(I)=MOD(IORIG(I),IPACK)    
        BP=BETA(1)*PPTCL(1,I)+BETA(2)*PPTCL(2,I)+BETA(3)*PPTCL(3,I) 
        DO 621 J=1,3    
          PPTCL(J,I)=PPTCL(J,I)+GAMMA*BETA(J)*(PPTCL(4,I)   
     $    +BP*GAMMA/(GAMMA+1.)) 
621     CONTINUE    
        PPTCL(4,I)=GAMMA*(PPTCL(4,I)+BP)    
620   CONTINUE  
C   
C          Normal exit  
C   
      DECJET=.TRUE. 
      RETURN    
C   
C          Error messages.  
C   
9998  DECJET=.FALSE.
      CALL PRTEVT(0)    
      WRITE(ITLIS,99980) NJSET  
99980 FORMAT(//5X,'ERROR IN DECJET...NJSET > ',I5)  
      RETURN    
      END   

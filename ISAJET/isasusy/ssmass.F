#include "isajet/pilot.h"
      SUBROUTINE SSMASS(XM1,XM2,IALLOW,ILOOP,MHLNEG,MHCNEG,IMODEL)
C-----------------------------------------------------------------------
C
C          Diagonalize neutralino, chargino, and Higgs mass matrices
C          and save results in /SSPAR/.
C
C          If XM1, XM2 < 1E19, use them for the U(1) and SU(2) mass
C          terms. Otherwise calculate them from AMGLSS and unification.
C
C          Return IALLOW = 1 if Z1SS is not LSP
C                 IALLOW = 0 otherwise
C
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
#include "isajet/sssm.inc"
#include "isajet/sspar.inc"
#include "isajet/ssinf.inc"
C
      REAL XM1,XM2
      INTEGER IALLOW,MHLNEG,MHCNEG,IMODEL
      REAL AR(4,4),WORK(4),WR(4)
      REAL ZETA,ZETAS,YM,XM,COS2A,SINA,AL,SIN2A,COSA,MU2,GP,G,
     $TEMP,VS,VP,V,MTAMTA,MTAMB,MTAMZ,ASMB,MBMB,
     $ASMT,MTMT,SUALFE,SUALFS
      REAL MW1,MW2,THX,THY,MU1
      REAL COSB,SINB,BE,COS2B,SIN2B,PI,SR2,HIGFRZ,MTQ,MBQ
      REAL TERM1,TERM2,TERM3,TANTHT,AMGLMZ,SSPOLE,TANTHB,TANTHL
      REAL CS2THW,DELCHI,AM2
      DOUBLE PRECISION SSMQCD
      COMPLEX*16 SSB0,SSB1,ZZZ
      REAL*8 REAL8
      INTEGER I,J,K,IERR,ILOOP
C
      REAL8(ZZZ)=DREAL(ZZZ)
      IALLOW=0
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4.*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      CS2THW=1.-SN2THW
C
      BE=ATAN(1./RV2V1)
      SINB=SIN(BE)
      COSB=COS(BE)
      SIN2B=SIN(2.*BE)
      COS2B=COS(2.*BE)
      HIGFRZ=MAX(AMZ,SQRT(AMTLSS*AMTRSS))
C
C          Compute m(tau), m(b) at z scale using qcd, qed
C
      MTAMTA=AMTAU*(1.-SUALFE(AMTAU**2)/PI)
      MTAMB=MTAMTA*(SUALFE(AMBT**2)/SUALFE(AMTAU**2))**(-27./76.)
      MTAMZ=MTAMB*(SUALFE(AMZ**2)/SUALFE(AMBT**2))**(-27./80.)
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SSMQCD(DBLE(MBMB),DBLE(HIGFRZ))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(HIGFRZ))
C
C     Light/heavy stop states and mixing angle
C
      TERM1=(AMTLSS**2+AMTRSS**2)/2.+AMZ**2*COS2B/4.+MTQ**2
      TERM2=((AMTLSS**2-AMTRSS**2)/2.+COS2B*(8.*AMW**2-5.*AMZ**2)
     $/12.)**2
      TERM3=SQRT(TERM2+MTQ**2*(TWOM1*COSB/SINB+AAT)**2)
      IF (TERM1.GT.TERM3) THEN
        AMT1SS=SQRT(TERM1-TERM3)
      ELSE
        AMT1SS=0.1
      END IF
      AMT2SS=SQRT(TERM1+TERM3)
      IF (AAT.NE.TWOM1*COSB/SINB) THEN
        TANTHT=(AMT1SS**2-MTQ**2+AMZ**2*COS2B*(-.5+2*SN2THW/3.)-
     $  AMTLSS**2)/MTQ/(TWOM1*COSB/SINB+AAT)
        THETAT=ATAN(TANTHT)
      ELSE
        THETAT=PI/2.
      END IF
C
C     Light/heavy sbottom states and mixing angle
C
      TERM1=(AMBLSS**2+AMBRSS**2)/2.-AMZ**2*COS2B/4.+MBQ**2
      TERM2=((AMBLSS**2-AMBRSS**2)/2.-COS2B*(4.*AMW**2-AMZ**2)
     $/12.)**2
      TERM3=SQRT(TERM2+MBQ**2*(TWOM1*SINB/COSB+AAB)**2)
      IF (TERM1.GT.TERM3) THEN
        AMB1SS=SQRT(TERM1-TERM3)
      ELSE
        AMB1SS=0.1
      END IF
      AMB2SS=SQRT(TERM1+TERM3)
      TANTHB=(AMB1SS**2-MBQ**2+AMZ**2*COS2B*(.5-SN2THW/3.)-
     $AMBLSS**2)/MBQ/(TWOM1*SINB/COSB+AAB)
      THETAB=ATAN(TANTHB)
C
C     Light/heavy stau states and mixing angle
C
      TERM1=(AMLLSS**2+AMLRSS**2)/2.-AMZ**2*COS2B/4.+MTAMZ**2
      TERM2=((AMLLSS**2-AMLRSS**2)/2.-COS2B*(4.*AMW**2-3*AMZ**2)
     $/4.)**2
      TERM3=SQRT(TERM2+MTAMZ**2*(TWOM1*SINB/COSB+AAL)**2)
C     if stau mass^2<0, then set to tiny mass so point is excluded
      IF (TERM1.GT.TERM3) THEN
        AML1SS=SQRT(TERM1-TERM3)
      ELSE
        AML1SS=0.1
      END IF
      AML2SS=SQRT(TERM1+TERM3)
      TANTHL=(AML1SS**2-MTAMZ**2+AMZ**2*COS2B*(.5-SN2THW)-
     $AMLLSS**2)/MTAMZ/(TWOM1*SINB/COSB+AAL)
      THETAL=ATAN(TANTHL)
C
C     define msbar gluino mass at mz from physical gluino mass
      AMGLMZ=SSPOLE(AMGLSS,AMZ**2,-ALFA3)
      VS=2.*AMW**2/G**2/(1.+RV2V1**2)
      V=SQRT(VS)
      VP=RV2V1*V
C
C          Use either explicit values or scaling to determine SU(2)
C          and U(1) mass terms. NOTE SIGN CONVENTION!
C
      IF(ABS(XM2).LT.1.E19.AND.ABS(XM1).LT.1.E19) THEN
         MU2=-XM2
         MU1=-XM1
      ELSE
         MU2=-ALFA2*AMGLMZ/ALFA3
         MU1=5*SN2THW/3./(1.-SN2THW)*MU2
      ENDIF
C
C          Neutralino mass matrix
C
      AR(1,1)=0.
      AR(1,2)=-TWOM1
      AR(1,3)=-G*V/SR2
      AR(1,4)=GP*V/SR2
      AR(2,1)=-TWOM1
      AR(2,2)=0.
      AR(2,3)=G*VP/SR2
      AR(2,4)=-GP*VP/SR2
      AR(3,1)=-G*V/SR2
      AR(3,2)=G*VP/SR2
      AR(3,3)=MU2
      AR(3,4)=0.
      AR(4,1)=GP*V/SR2
      AR(4,2)=-GP*VP/SR2
      AR(4,3)=0.
      AR(4,4)=MU1
C
      CALL EISRS1(4,4,AR,WR,ZMIXSS,IERR,WORK)
      IF (IERR.NE.0) THEN
        WRITE(LOUT,*) 'EISRS1 ERROR IN SSMASS, IERR=',IERR
        STOP99
      END IF
C
C       Sort eigenvectors and eigenvalues according to masses
C
      DO 10 I=1,3
        DO 11 J=I+1,4
          IF (ABS(WR(I)).GT.ABS(WR(J))) THEN
            TEMP=WR(J)
            WR(J)=WR(I)
            WR(I)=TEMP
            DO 12 K=1,4
              TEMP=ZMIXSS(K,J)
              ZMIXSS(K,J)=ZMIXSS(K,I)
              ZMIXSS(K,I)=TEMP
12          CONTINUE
          END IF
11      CONTINUE
10    CONTINUE
C
      AMZ1SS=WR(1)
      AMZ2SS=WR(2)
      AMZ3SS=WR(3)
      AMZ4SS=WR(4)
C
C          Chargino mass matrix
C
      AL=ATAN(RV2V1)
      SINA=SIN(AL)
      COSA=COS(AL)
      SIN2A=SIN(2.*AL)
      COS2A=COS(2.*AL)
      ZETAS=(TWOM1**2-MU2**2)**2
     $+4*AMW**2*(AMW**2*COS2A**2+TWOM1**2+MU2**2+2*TWOM1*MU2*SIN2A)
      ZETA=SQRT(ZETAS)
      XM=-(TWOM1**2-MU2**2-2*AMW**2*COS2A-ZETA)
     $/(2*SR2*AMW*(MU2*SINA+TWOM1*COSA))
      YM=-(TWOM1**2-MU2**2+2*AMW**2*COS2A-ZETA)
     $/(2*SR2*AMW*(MU2*COSA+TWOM1*SINA))
      IF (XM.NE.0.) THEN
        GAMMAL=ATAN(1./XM)
      ELSE
        GAMMAL=PI/2.
      END IF
      IF (YM.NE.0.) THEN
        GAMMAR=ATAN(1./YM)
      ELSE
        GAMMAR=PI/2.
      END IF
      IF (GAMMAL.LT.0.) GAMMAL=GAMMAL+PI
      IF (GAMMAR.LT.0.) GAMMAR=GAMMAR+PI
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      AMW2SS=THX*THY*(COS(GAMMAR)*(MU2*COS(GAMMAL)+G*VP*SIN(GAMMAL))
     $-SIN(GAMMAR)*(-G*V*COS(GAMMAL)-TWOM1*SIN(GAMMAL)))
      AMW1SS=SIN(GAMMAR)*(MU2*SIN(GAMMAL)-G*VP*COS(GAMMAL))
     $+COS(GAMMAR)*(-G*V*SIN(GAMMAL)+TWOM1*COS(GAMMAL))
C     IMPLEMENT INO MASS SPLITTING FOR AMSB MODELS
      AM2=ABS(XM2)
      XLAM=LOG(MU2**2)
      MW1=ABS(AMW1SS)
      DELCHI=G**2*MW1/8./PI**2*(2*CS2THW*REAL8(SSB0(MW1**2,MW1,AMZ))+
     $2*SN2THW*REAL8(SSB0(MW1**2,MW1,0.))-2*REAL8(SSB0(MW1**2,MW1,AMW))
     $-CS2THW*REAL8(SSB1(MW1**2,MW1,AMZ))-SN2THW*
     $REAL8(SSB1(MW1**2,MW1,0.))+REAL8(SSB1(MW1**2,MW1,AMW)))
      AMW1SS=AMW1SS+SIGN(1.,AMW1SS)*DELCHI
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
C
C          Check validity of parameters
C
      IF (IMODEL.EQ.1.OR.IMODEL.EQ.7) THEN
        IF(MW1.LE.ABS(AMZ1SS)) IALLOW=1
        IF(AMT1SS.LE.ABS(AMZ1SS)) IALLOW=1
        IF(AMB1SS.LE.ABS(AMZ1SS)) IALLOW=1
        IF(AML1SS.LE.ABS(AMZ1SS)) IALLOW=1
      END IF
C      IF(IALLOW.NE.0) RETURN
C
C          Higgs mass matrix
C
      IF (ILOOP.EQ.1) THEN
        CALL SSMHN(MHLNEG)
        CALL SSMHC(MHCNEG)
      END IF
C
      RETURN
      END

      SUBROUTINE CPVREC(DzCPV,DyCPV,DxCPV,Nhit,Xin,Pin,Nout,Xout)
C-----------------------------------------------------------------------
C     Reconstruction routing for the CPV detector
C
C     Author: Serguei Sadovski, IHEP, Protvino
C     Created:       25 December 1998
C     Last modified: 17 September 1999
C-----------------------------------------------------------------------
C
C CPV Frame:  Z coordinate - along beam direction
C ==========  Y coordinate - transverse beam direction in the CPV plane
C             X coordinate - transverse CPV plane, positive direction - along
C                            particles going from the interaction point
C            (0,0) in the center of the CPV plane
C
C Input Parameters:
C =================
C DzCPV     - half width (cm) of CPV plane along beam
C DyCPV     - half width (cm) of CPV plane transverse beam
C DxCPV     - half width (cm) of CPV
C Nhit      - number of particles in the front of CPV
C
C Xin(1,i)  - Z coordinate (cm) of i-hit in CPV frame (0,0 in the CPV Center)
C Xin(2,i)  - Y coordinate (cm) of i-hit in CPV frame (0,0 in the CPV Center)
C
C Pin(1,i)  - Pz
C Pin(2,i)  - Py      for particle i
C Pin(3,i)  - Px
C Pin(4,i)  - dE/dX
C
C Output Parameters - reconstructed coordinates of particles in the CPV:
C =================
C Nout      - number of the reconstructed hits (<5000) in CPV plane,
C             if (Nout.le.0) - error output: reconstruction fails ...
C
C Xout(1,j) - Z coordinate (cm) of j-hit in CPV frame (0,0 in the CPV Center)
C Xout(2,j) - Y coordinate (cm) of j-hit in CPV frame (0,0 in the CPV Center)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DIMENSION Xin(2,Nhit),Pin(4,Nhit),Xout(2,5000)             ! Inp/Out
C
      COMMON/GAMMA/ NGAM,EGAM(5000),XGAM(5000),YGAM(5000),CHGAM(5000)! From Rec
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
C
      REAL   *4 EGM,ZGM,YGM, ZYEgm(3,5)
      REAL   *8 Pp
      INTEGER*2               ICPV(104,176)                      ! CPV 104x176
      DATA   NcelZ,NcelY,CelZ,CelY/104,176, 2.26, 1.13 /         ! in cm
      DATA   NgamZ,NgamY,Chit,Cnoise/5,9  , 1000., 30. /
      DATA   CelWr,RhoWr,IDELT,IDELA/ 0.5   , .01, 90,  90 /     ! 09.01.99
      DATA   Detr / 0.1 /  !  Relative energy fluctuation in Track for 100 e-
C
      Nout = 0
      Ncpv = 0
      IF (Nhit.eq.0) RETURN
C!-                              - - -  C O N S T A N T S  - - -
      ECPV = 0.
      NZMX = 0
      NZMN = NcelZ
      NYMX = 0
      NYMN = NcelY
C
      CALL VZERO(ICPV,NcelZ*NcelY/2)
C
      DO 1000 Ihit=1,Nhit
C
      Pp   =  SQRT( Pin(1,Ihit)**2 +  Pin(2,Ihit)**2 +  Pin(3,Ihit)**2 )
C
      If (Abs(Pin(3,Ihit))/Pp.lt.0.0000001) go to 1000 ! Tracks paralell CPV
C
      CorX=    2.*DxCPV/Pin(3,Ihit)
      Dzx = Pin(1,Ihit)*CorX            !  Incline tracks
      Dyx = Pin(2,Ihit)*CorX
C                                       !  Charge of inline tracks
      CALL RNORML(R1,1)
      Pin(4,Ihit) =  Pin(4,Ihit)*(1.+Detr*R1)*
     +               SQRT(Dzx**2 + Dyx**2 + 4.*DxCPV**2)/(2.*DxCPV)
C
c      call hfill(113,Dzx, 0., 1.)
c      call hfill(114,Dyx, 0., 1.)
C
      Zhit1 = Xin(1,Ihit)+DzCPV
      Yhit1 = Xin(2,Ihit)+DyCPV         !  Initial Y coordinate
C
      Zhit2 = Zhit1+Dzx                 !  final Y coordinate
      Yhit2 = Yhit1+Dyx
C
      Iwht1 = Yhit1/CelWr               !  Wire (Y) coordinate "in"
      Iwht2 = Yhit2/CelWr               !  Wire (Y) coordinate "out"
C
      If(Iwht1.eq.Iwht2)   then         !  Incline 1 wire hit
                           Niter      = 2
                           ZYEgm(1,1) =(Zhit1+Zhit2-Dzx*.57735)/2.
                           ZYEgm(2,1) =(Iwht1+0.5)*CelWr           !
                           ZYEgm(3,1) = Pin(4,Ihit)/2.
C
                           ZYEgm(1,2) =(Zhit1+Zhit2+Dzx*.57735)/2.
                           ZYEgm(2,2) =(Iwht1+0.5)*CelWr           !
                           ZYEgm(3,2) = Pin(4,Ihit)/2.
                        go to 10
                  endif
C                                       !    3 Wire tracks   !
C                                     --------------------------
      IF(Iwht1-1.ne.Iwht2 .and.
     +   Iwht1+1.ne.Iwht2)     THEN
                               Niter = 3
                               Iwht3 =(Iwht1+Iwht2)/2
                               EGM   = Pin(4,Ihit)
                               Ywht1 =(Iwht1+0.5)*CelWr          !  Wire 1
                               Ywht2 =(Iwht2+0.5)*CelWr          !  Wire 2
                               Ywht3 =(Iwht3+0.5)*CelWr          !  Wire 3
                               Ywr13 =(Ywht1+Ywht3)/2.           !  Center 13
                               Ywr23 =(Ywht2+Ywht3)/2.           !  Center 23
C
           DyW1  = Yhit1-Ywr13
           DyW2  = Yhit2-Ywr23
           Egm1  = ABS(DyW1)/(ABS(DyW1)+ABS(DyW2) + CelWr)
           Egm2  = ABS(DyW2)/(ABS(DyW1)+ABS(DyW2) + CelWr)
           Egm3  =    CelWr /(ABS(DyW1)+ABS(DyW2) + CelWr)
C
           ZYEgm(1,1) =(Dzx*(Ywr13-Ywht1)/Dyx + Zhit1 + Zhit1)/2.
           ZYEgm(2,1) = Ywht1
           ZYEgm(3,1) = EGM*Egm1
C
           ZYEgm(1,2) =(Dzx*(Ywr23-Ywht1)/Dyx + Zhit1 + Zhit2)/2.
           ZYEgm(2,2) = Ywht2
           ZYEgm(3,2) = EGM*Egm2
C
           ZYEgm(1,3) = Dzx*(Ywht3-Ywht1)/Dyx + Zhit1
           ZYEgm(2,3) = Ywht3
           ZYEgm(3,3) = EGM*Egm3
C
           go to 10
      ENDIF
C
      EGM   = Pin(4,Ihit)
      Ywht1 =(Iwht1+0.5)*CelWr          !  Wire 1
      Ywht2 =(Iwht2+0.5)*CelWr          !  Wire 2
      Ywr12 =(Ywht1+Ywht2)/2.
C
      DyW1  = Yhit1-Ywr12
      DyW2  = Yhit2-Ywr12
      Egm1  = ABS(DyW1)/(ABS(DyW1)+ABS(DyW2))
      Egm2  = ABS(DyW2)/(ABS(DyW1)+ABS(DyW2))
C
      Niter      = 2
      ZYEgm(1,1) =(Zhit1+Zhit2-Dzx*Egm1)/2.
      ZYEgm(2,1) = Ywht1
      ZYEgm(3,1) = EGM*Egm1
C
      ZYEgm(1,2) =(Zhit1+Zhit2+Dzx*Egm2)/2.
      ZYEgm(2,2) = Ywht2
      ZYEgm(3,2) = EGM*Egm2
C
c      CALL HFILL(110, Egm1,0.,1.)
c      CALL HFILL(110, Egm2,0.,1.)
C
C                                       !  Finite size of ionisation region
  10  DO 113 Iter=1,Niter
      Zhit =  ZYEgm(1,Iter)             !
      Yhit =  ZYEgm(2,Iter)             !  Wire coordinate
      EGM  =  ZYEgm(3,Iter)
C
      IF   ( Zhit.LE.0. .OR. Yhit.LE.0. ) GO TO 1000
C
      Zcel = Zhit/CelZ
      Ycel = Yhit/CelY
      Izg  = Zcel
      Iyg  = Ycel
      IF ( Izg.GE.NcelZ.OR.Iyg.GE.NcelY ) GO TO 1000
C!-
      Zc   = Zcel - Izg-.5
      Yc   = Ycel - Iyg-.5
C
      Nz3  = (NgamZ+1)/2
      Ny3  = (NgamY+1)/2
      DO 112 JJ=1,NgamZ
      J=JJ-Nz3
      KZG=IZG+J
      IF (KZG.LE.0.OR.KZG.GT.NcelZ) GO TO 112
      ZG=FLOAT(J)-Zc
C
      DO 111 II=1,NgamY
      I=II-Ny3
      KYG=IYG+I
      IF (KYG.LE.0.OR.KYG.GT.NcelY) GO TO 111
      YG=FLOAT(I)-Yc
C
      CALL RNORML(R1,1)
      AEG = AMP(EGM,ZG,YG)
      AGM = Chit*AEG + R1*Cnoise
      IF  ( AGM.LE.0.)              GO TO 111
C
      ICPV(KZG,KYG) = ICPV(KZG,KYG) + AGM
      ECPV  = ECPV  +                 AGM
C
C-    type *, ' KZG,KYG,AGM=',KZG,KYG,AGM
C
      IF  (KZG.LE.NZMN) NZMN = KZG
      IF  (KZG.GE.NZMX) NZMX = KZG
      IF  (KYG.LE.NYMN) NYMN = KYG
      IF  (KYG.GE.NYMX) NYMX = KYG
 111  CONTINUE
 112  CONTINUE
 113  CONTINUE
      Ncpv = Ncpv + 1
C
 1000 CONTINUE
c      CALL HF1( 2,ECPV,1.)
c      CALL HF1(31,ECPV,1.)
C
      CALL COMPRS(NZMN,NcelZ,NYMN,NcelY,ICPV,NWGAMS)
      CALL RECONS(NWGAMS)
C
      IF(NGAM.LE.0) RETURN
C                                                !   O
      DO Ig=1,NGAM                               !   U
      Xout(1,Ig)  = XGAM(Ig) + 0.5*CelZ - DzCPV  !   T
      Xout(2,Ig)  = YGAM(Ig) + 0.5*CelY - DyCPV  !   P
      ENDDO                                      !   U
      Nout = NGAM                                !   T
C
      RETURN
      END
C=======================================================================
      SUBROUTINE COMPRS(NXMN,NXMX,NYMN,NYMX,IA,NWGAMS)
      COMMON/EVENT/ IBH(5),IPHD,IPSC,IPGS,IPGM,IPTG,I11,I12,IBUF(30035)
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
      INTEGER*2  IPH,IA(NXMX,NYMX)
      DATA  IPGM,Mdim,IBUF(3) / 35, 30035, 0 /
C
      IPGMS =IPGM
      DO 10 IY=NYMN,NYMX
      DO 10 IX=NXMN,NXMX
C
      IPH = IA(IX,IY)
      IF (IPH.LT.IDELT) GO TO 10
c      CALL HFILL(9,FLOAT(IPH),0.,1.)
C
      IPGMS =IPGMS+2
      IF(IPGMS.ge.Mdim) GO TO 10
C
      ICNT  =NYMX* IX   +IY
      IBUF(  IPGMS) = IPH
      IBUF(1+IPGMS) = ICNT
   10 CONTINUE
C
      NWGAMS = IPGMS-IPGM

      IBUF(1)= IPGMS
      RETURN
      END
C
C=======================================================================
      SUBROUTINE RECONS(NWGAMS)
      COMMON/EVENT/ IHDR(12),IBUF(30035)
      EQUIVALENCE(IHDR(3),NREC),(IHDR(5),NMT),(IHDR(6),IPHOD),(IHDR(7),
     ,IPSC),(IHDR(8),IPGS),(IHDR(9),IPGM),(IHDR(10),IPTG),(IBUF(3),NEV)
      COMMON/GAMMA/ NGAM,EGAM(5000),XGAM(5000),YGAM(5000),CHGAM(5000)
      COMMON/CLUSTER/ NCL,LENCL(5000)
      COMMON/WORK/ IWRK(40000)
      DATA NWMIN,NWMAX/2,30000/,MINPK,MINECL/10,15/,ECUT/25./
      DATA SCALMS/26./
C
      NGAM=0
c      CALL HFILL(1,  FLOAT(NWGAMS),  0.,1.)   !..my..!
c      CALL HFILL(101,FLOAT(NWGAMS/2),0.,1.)

      IF(NWGAMS.LT.NWMIN.OR.NWGAMS.GT.NWMAX) RETURN
      IPGAMS=IPGM+2
      ETOT=0.
      IPEND=IPGAMS+NWGAMS-2
      DO 1 I=IPGAMS,IPEND,2
    1 ETOT=ETOT+IBUF(I)
c      CALL HFILL( 7,ETOT,0.,1.)
c      CALL HFILL(32,ETOT,0.,1.)
      IF(ETOT.LT.ECUT)         RETURN

      CALL CLUST(NWGAMS,IBUF(IPGAMS))
c      CALL HFILL(3,FLOAT(NCL),0.,1.)  !..my..!
C
      CALL VARDEL(NWGAMS)
      IF(NCL.LT.1)             RETURN
      ETOTCL=0.
      IPCL=IPGAMS
      DO 10 ICL=1,NCL
         IECL=0
C     
c         CALL HFILL(4,FLOAT(LENCL(ICL)),0.,1.) !..my..!
C     
         DO 5 I=1,LENCL(ICL)
 5       IECL=IECL+IBUF(IPCL+2*I-2)
         IF (IECL.GE.MINECL)       ETOTCL=ETOTCL+IECL
*
c         IF (iecl.GE.900) THEN
c            CALL hf1 (104,float(lencl(icl)),1.)
c         END IF
*
         IPCL=IPCL+LENCL(ICL)*2
 10   CONTINUE
c      CALL HFILL( 5,ETOTCL,0.,1.)  !..my..!
c      CALL HFILL(33,ETOTCL,0.,1.)  !..my..!
      IF(ETOTCL.LT.ECUT)       RETURN
C
      IPCL=IPGAMS
      DO 30 ICL=1,NCL
      IECL=0
      DO 25 I=1,LENCL(ICL)
   25 IECL=IECL+IBUF(IPCL+2*I-2)
      IF(IECL.LT.MINECL)       GO TO 30
      CALL GAMS(IBUF(IPCL),LENCL(ICL),IWRK(1),MINPK)
   29 IF(NGAM.GE.4999)         GO TO 50
   30 IPCL=IPCL+LENCL(ICL)*2
C
   50 CONTINUE
      CALL KILLG(SCALMS)
c      CALL HFILL(6,FLOAT(NGAM),0.,1.)  !..my..!
      RETURN
      END
C=======================================================================
      SUBROUTINE KILLG(SCALMS)
      COMMON /GAMMA / NGAM,EGAM(5000),XGAM(5000),YGAM(5000),CHGAM(5000)
      COMMON /CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelYC
cc    DATA   SQDCUT / 3.50/, EFMCUT / 5./, THR0 /500./  ! N_c = 1500
      DATA   SQDCUT / 4.50/, EFMCUT / 5./, THR0 /600./  ! N-c = 800
C
      NGM=0
      ETOT=0.
      DO 1 IG=1,NGAM
         IF(EGAM(IG).LT.THR0)                 GO TO 1
         NGM=NGM+1
         XGAM(NGM) = XGAM(IG)
         YGAM(NGM) = YGAM(IG)
         EGAM(NGM) = EGAM(IG)
         CHGAM(NGM)=CHGAM(IG)
         ETOT=ETOT + EGAM(IG)
c         CALL HF1(35,EGAM(NGM),1.)
    1 CONTINUE
      NGAM=NGM
c      CALL HF1(37,FLOAT(NGM),1.)
C
      iloop = 0
    2 IF(NGAM.LT.2)                           RETURN
      SQDMIN = 1.0E+10
      iloop = iloop + 1
      DO IG1=1,NGAM-1
         DO IG2=IG1+1,NGAM
            SQDIST = (XGAM(IG1)-XGAM(IG2))**2 + (YGAM(IG1)-YGAM(IG2))**2
c            IF (iloop .EQ. 1) THEN
c               CALL HF1(36, SQDIST, 1.)
c            END IF
            IF (SQDIST .LT. SQDMIN) THEN
               SQDMIN = SQDIST
               IG1M   = IG1
               IG2M   = IG2
            END IF
         END DO
      END DO
C
      IF (SQDMIN .GT. SQDCUT) RETURN
*
      IF (egam(ig1m) .GT. egam(ig2m)) THEN
         XGAM(IG1M)  =  XGAM(IG1M)
         YGAM(IG1M)  =  YGAM(IG1M)
         CHGAM(IG1M) = CHGAM(IG1M)
         EGAM(IG1M)  =  EGAM(IG1M)
      ELSE
         XGAM(IG1M)  =  XGAM(IG2M)
         YGAM(IG1M)  =  YGAM(IG2M)
         CHGAM(IG1M) = CHGAM(IG2M)
         EGAM(IG1M)  =  EGAM(IG2M)
      END IF
c      DE          =  EGAM(IG1M)/(EGAM(IG1M)+EGAM(IG2M))
c      XGAM(IG1M)  =  XGAM(IG1M)*DE +  XGAM(IG2M)*(1.-DE)
c      YGAM(IG1M)  =  YGAM(IG1M)*DE +  YGAM(IG2M)*(1.-DE)
c      CHGAM(IG1M) = CHGAM(IG1M)    + CHGAM(IG2M)
c      EGAM(IG1M)  =  EGAM(IG1M)    +  EGAM(IG2M)

c      CALL hf1 (38, egam(ig1m), 1.)
      NGAM=NGAM-1
      IF(IG2M.GT.NGAM)                        GO TO 2
      DO IG=IG2M,NGAM
         XGAM(IG) = XGAM(IG+1)
         YGAM(IG) = YGAM(IG+1)
         EGAM(IG) = EGAM(IG+1)
         CHGAM(IG)=CHGAM(IG+1)
      END DO
      GO TO 2
*
      END
*=======================================================================
      SUBROUTINE ORDER(N,IAR)
      DIMENSION IAR(N)
      IF(N.LT.4) RETURN
      DO 4 K=4,N,2
      IF(IAR(K).GE.IAR(K-2)) GO TO 4
      IA=IAR(K)
      ID=IAR(K-1)
      IF(K.EQ.4) GO TO 2
      DO 1 I1=2,K-4,2
      I=K-I1
      IF(IAR(I-2).LE.IA) GO TO 3
      IAR(I+2)=IAR(I)
    1 IAR(I+1)=IAR(I-1)
    2 I=2
    3 IAR(I+2)=IAR(I)
      IAR(I+1)=IAR(I-1)
      IAR(I)=IA
      IAR(I-1)=ID
    4 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE CLUST(NW,IA)
      DIMENSION IA(2,NW)
      COMMON/CLUSTER/ NCL,LENCL(5000)
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
      COMMON/WORK/ IWORK(40000)
      DATA MxCL  / 5000 /
      NCL=0
      IF(NW.LT.2) RETURN
      NCL=1
      LENCL(1)=1
      IF(NW.LT.4) RETURN
      CALL ORDER(NW,IA(1,1))
      NW2=NW/2
      NCL=0
      NEXT=IA(2,1)
      IAK=NEXT
      K2 =2
    5 DO 10 K=K2,NW2
      IAK1=IAK
      IAK =IA(2,K)
C
      IF(IAK.EQ.IAK/NcelY*NcelY) GO TO 20
      IF(IAK.GT.IAK1+1)          GO TO 20
   10 CONTINUE
C
   20 IBEG=NEXT
      NEXT=IAK
      IEND=IA(2,K-1)
      KB  = K-1 - IEND + IBEG
      IF(NCL.GE.MxCL)            RETURN
C
      NCL=NCL+1
      LENCL(NCL)=IEND-IBEG+1
      IF(NCL.EQ.1.AND.K.GT.NW2)  RETURN
      IF(NCL.EQ.1)   THEN
                     K2=K+1
                     GO TO 5
                     ENDIF
      LAST=KB-1
      NCL0=NCL
      DO 60 ICL1=1,NCL0-1
      ICL=NCL0-ICL1
      IF(IBEG-IA(2,LAST).GT.NcelY.AND.K.GT.NW2) RETURN
      IF(IBEG-IA(2,LAST).GT.NcelY)  THEN
                                    K2=K+1
                                    GO TO 5
                                    ENDIF
      LENG=LENCL(ICL)
C-
      DO 50 I=1,LENG
      ICUR=IA(2,LAST-I+1)
      IF(IEND-ICUR.LT.NcelY)        GO TO 50
      IF(IBEG-ICUR.GT.NcelY)        GO TO 60
      IF(ICL.EQ.NCL-1)           GO TO 40
      IF(LENG.GT.768)            GO TO 60
      CALL UCOPY(IA(1,LAST+1-LENG),IWORK(1),LENG*2)
      CALL UCOPY(IA(1,LAST+1),IA(1,LAST+1-LENG),(KB-1-LAST)*2)
      CALL UCOPY(IWORK(1),IA(1,KB-LENG),LENG*2)
      DO 30 J=ICL,NCL-2
   30 LENCL(J)=LENCL(J+1)
   40 NCL=NCL-1
      KB = KB - LENG
      LENCL(NCL)=K-KB
      GO TO 60
   50 CONTINUE
   60 LAST=LAST-LENG
      IF(K.LE.NW2)   THEN
                     K2=K+1
                     GO TO 5
                     ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE GAMS(IGAMS,NADC,IWRK,MINPK)
C  15.MAR.1983.
      COMMON/GAMMA/ NGAM,EGAM(5000),XGAM(5000),YGAM(5000),CHGAM(5000)
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
      DIMENSION IPNPK(50),EPK(50),XPK(50),YPK(50),IGMPK(2,50)
      DIMENSION IGAMS(2,NADC),IWRK(NADC,50)
      DATA MxPK                       / 50/,CHISQ1/30./,CHISQ2/3./
C
      IDELTA = IDELA
C
      CALL ORDER(NADC*2,IGAMS(1,1))
      NGAM0=NGAM
C
C   PEAKS SEARCH
C
      NPK=0
      DO 50 IC=1,NADC
      IAC=IGAMS(1,IC)
      IF(IAC.LT.MINPK)       GO TO 50
      IXY=IGAMS(2,IC)
      IXYC=IXY
      IYC=IXYC-IXYC/NcelY*NcelY
      IF(NADC.LT.2)          GO TO 40
      IF(IC.EQ.NADC)         GO TO 20
      DO 10 IN=IC+1,NADC
      IXY=IGAMS(2,IN)
      IF(IXY-IXYC.GT.NcelY+1)    GO TO 20
      IF(IABS(IXY-IXY/NcelY*NcelY-IYC).GT.1)  GO TO 10
      IF(IGAMS(1,IN).GE.IAC) GO TO 50
   10 CONTINUE
   20 IF(IC.LT.2)            GO TO 40
      DO 30 IN1=1,IC-1
      IN=IC-IN1
      IXY=IGAMS(2,IN)
      IF(IXYC-IXY.GT.NcelY+1)     GO TO 40
      IF(IABS(IXY-IXY/NcelY*NcelY-IYC).GT.1)  GO TO 30
      IF(IGAMS(1,IN).GT.IAC) GO TO 50
   30 CONTINUE
   40 NPK=NPK+1
      IPNPK(NPK)=IC
C!!-  IF(NPK.EQ.  10.OR.NPK.GE.NcelY*NcelY/NADC-3) GO TO 60
      IF(NPK.EQ.MxPK.OR.NPK.GE.NcelZ*NcelY/NADC-3) GO TO 60
   50 CONTINUE
      IF(NPK.EQ.0)           RETURN
C
   60 CONTINUE
C
C   GAMMA SEARCH FOR ONE PEAK
C
  101 IF(NPK.GT.1)           GO TO 102
      IF(NGAM.GE.4999)       RETURN
      NGAM=NGAM+1
      CHISQ=CHISQ2
      CALL GAM (NADC,IGAMS(1,1),CHISQ,EGAM(NGAM),XGAM(NGAM),YGAM(NGAM),
     ,                          EGAM(NGAM+1),XGAM(NGAM+1),YGAM(NGAM+1))
      CHGAM(NGAM)=CHISQ
      IF(EGAM(NGAM+1).EQ.0..OR.NGAM.EQ.4999)  GO TO 500
      NGAM=NGAM+1
      CHGAM(NGAM)=CHISQ
      GO TO 500
C
C   GAMMA SEARCH FOR SEVERAL PEAKS
C   FIRST STEP AND STEP 1-BIS.
C
  102 IF(NGAM.GE.4999)         RETURN
      DO 170 IBIS=1,2
      DO 105 I=1,NADC
  105 IWRK(I,NPK+3)=0
C
      DO 160 IPK=1,NPK
      IC=IPNPK(IPK)
      E=IGAMS(1,IC)
      IF(IBIS.NE.1) E=E*FLOAT(IWRK(IC,IPK))/IWRK(IC,NPK+1)
      IXYPK=IGAMS(2,IC)
      IXPK=IXYPK/NcelY
      IYPK=IXYPK-IXPK*NcelY
      EPK(IPK)=E
      XPK(IPK)=E*IXPK
      YPK(IPK)=E*IYPK
      IF(IC.GE.NADC)         GO TO 120
      DO 110 IN=IC+1,NADC
      IXY=IGAMS(2,IN)
      IF(IXY-IXYPK.GT.NcelY+1)    GO TO 120
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      IF(IABS(IY-IYPK).GT.1) GO TO 110
      E=IGAMS(1,IN)
      IF(IBIS.NE.1) E=E*FLOAT(IWRK(IN,IPK))/IWRK(IN,NPK+1)
      EPK(IPK)=EPK(IPK)+E
      XPK(IPK)=XPK(IPK)+E*IX
      YPK(IPK)=YPK(IPK)+E*IY
  110 CONTINUE
  120 IF(IC.LT.2)            GO TO 140
      DO 130 IN1=1,IC-1
      IN=IC-IN1
      IXY=IGAMS(2,IN)
      IF(IXYPK-IXY.GT.NcelY+1)    GO TO 140
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      IF(IABS(IY-IYPK).GT.1) GO TO 130
      E=IGAMS(1,IN)
      IF(IBIS.NE.1) E=E*FLOAT(IWRK(IN,IPK))/IWRK(IN,NPK+1)
      EPK(IPK)=EPK(IPK)+E
      XPK(IPK)=XPK(IPK)+E*IX
      YPK(IPK)=YPK(IPK)+E*IY
  130 CONTINUE
  140 CONTINUE
      XPK(IPK)=XY(XPK(IPK)/EPK(IPK))
      YPK(IPK)=XY(YPK(IPK)/EPK(IPK))
      DO 150 I=1,NADC
      IXY=IGAMS(2,I)
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      DX=IX-XPK(IPK)
      DY=IY-YPK(IPK)
      IWK=0
      IF(DX*DX+DY*DY.LT.8.) IWK = AMP(EPK(IPK),DX,DY)+.5
      IWRK(I,IPK)=IWK
      IWRK(I,NPK+3)=IWRK(I,NPK+3)+IWK
  150 CONTINUE
  160 CONTINUE
      DO 165 I=1,NADC
      IWK=IWRK(I,NPK+3)
      IF(IWK.LE.0) IWK=1
  165 IWRK(I,NPK+1)=IWK
  170 CONTINUE
C
      DO 200 IPK=1,NPK
      LENG=0
      DO 190 I=1,NADC
      IF(IWRK(I,NPK+3).LE.0) GO TO 190
      IE=IGAMS(1,I)*FLOAT(IWRK(I,IPK))/IWRK(I,NPK+3)+.5
      IF(IE.LE.IDELTA)       GO TO 190
      LENG=LENG+1
      IWRK(2*LENG-1,NPK+1)=IE
      IWRK(2*LENG,NPK+1)=IGAMS(2,I)
  190 CONTINUE
      IF(NGAM.GE.4999)       GO TO 500
      IF(LENG.EQ.0)          GO TO 200
      NGAM=NGAM+1
      CHISQ=CHISQ1
      CALL GAM (LENG,IWRK(1,NPK+1),CHISQ,EGAM(NGAM),XGAM(NGAM),
     ,       YGAM(NGAM),EGAM(NGAM+1),XGAM(NGAM+1),YGAM(NGAM+1))
      CHGAM(NGAM)=CHISQ
      IGMPK(1,IPK)=NGAM
      IGMPK(2,IPK)=NGAM
      IF(EGAM(NGAM+1).EQ.0..OR.NGAM.EQ.4999)   GO TO 200
      NGAM=NGAM+1
      CHGAM(NGAM)=CHISQ
      IGMPK(2,IPK)=NGAM
  200 CONTINUE
C
C   SECOND STEP.
C
      DO 210 I=1,NADC
  210 IWRK(I,NPK+3)=0
      DO 230 IPK=1,NPK
      DO 230 I=1,NADC
      IWRK(I,IPK)=0
      DO 220 IG=IGMPK(1,IPK),IGMPK(2,IPK)
      IXY=IGAMS(2,I)
      DX=IXY/NcelY-XGAM(IG)
      DY=IXY-IXY/NcelY*NcelY-YGAM(IG)
      IF(DX*DX+DY*DY.GT.8.)  GO TO 220
      IWK = AMP(EGAM(IG),DX,DY)+.5
      IWRK(I,IPK)  =IWRK(I,IPK)  +IWK
      IWRK(I,NPK+3)=IWRK(I,NPK+3)+IWK
  220 CONTINUE
  230 CONTINUE
C
      NGAM=NGAM0
      DO 300 IPK=1,NPK
      LENG=0
      DO 290 I=1,NADC
      IF(IWRK(I,NPK+3).LE.0) GO TO 290
      IE=IGAMS(1,I)*FLOAT(IWRK(I,IPK))/IWRK(I,NPK+3)+.5
      IF(IE.LE.IDELTA)       GO TO 290
      LENG=LENG+1
      IWRK(2*LENG-1,NPK+1)=IE
      IWRK(2*LENG,NPK+1)=IGAMS(2,I)
  290 CONTINUE
      IF(NGAM.GE.4999)       GO TO 500
      IF(LENG.EQ.0)          GO TO 300
      NGAM=NGAM+1
      CHISQ=CHISQ2
      CALL GAM (LENG,IWRK(1,NPK+1),CHISQ,EGAM(NGAM),XGAM(NGAM),
     ,       YGAM(NGAM),EGAM(NGAM+1),XGAM(NGAM+1),YGAM(NGAM+1))
      CHGAM(NGAM)=CHISQ
      IF(EGAM(NGAM+1).EQ.0..OR.NGAM.EQ.4999)   GO TO 300
      NGAM=NGAM+1
      CHGAM(NGAM)=CHISQ
  300 CONTINUE
  500 DO 600 IG=NGAM0+1,NGAM
      XGAM(IG)= XGAM(IG)*CelZ
  600 YGAM(IG)= YGAM(IG)*CelY
      RETURN
      END
C=======================================================================
      SUBROUTINE GAM(NADC,IADC,CHISQ,E1,X1,Y1,E2,X2,Y2)
C   28 OCT.1981.     ONE GAMMA FIT PROGRAM.
CORRECTED 28 MAR.1983.
      DIMENSION IADC(2,NADC)
      DATA STPMIN/.005/,CONST/1.5/,ST0/.00005/,CHMIN/1./
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
C
C     WRITE(5,1)
C   1 FORMAT(' *** GAM ***')
      IF(NADC.LE.0)                   RETURN
      DOF=NADC-2
      IF(DOF.LT.1.)                   DOF=1.
      CHSTOP=CHMIN*DOF
      CSQCUT=CHISQ
      E1=0.
      XX=0.
      YY=0.
      E2=0.
      X2=0.
      Y2=0.
      DO 5 I=1,NADC
      A  = IADC(1,I)
      E1 = E1 + A
      IXY=IADC(2,I)
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      XX = XX + A*IX
      YY = YY + A*IY
    5 CONTINUE
      XC = XY(XX/E1)
      YC = XY(YY/E1)
      ST = ST0
      CHISQ=1.E20
      GR = 1.
      GRX=0.
      GRY=0.
      DO 100 ITER=1,50
      CHISQC=0.
      GRXC  =0.
      GRYC  =0.
      DO 10 I=1,NADC
      IXY=IADC(2,I)
      DX = XC - IXY/NcelY
      DY = YC - (IXY-IXY/NcelY*NcelY)
      CALL AG(E1,DX,DY,A,GX,GY)
      D  = CONST * A * (1.-A/E1)
      IF(D.LT.0.) D = 0.
      D  =     9.
      EX = IADC(1,I)
      DA = A - EX
      CHISQC=CHISQC+DA**2/D
      DD = DA/D
      DD = DD*(2.-DD*CONST*(1.-2.*A/E1))
      GRXC = GRXC + GX*DD
      GRYC = GRYC + GY*DD
   10 CONTINUE
      GRC  = SQRT(GRXC**2+GRYC**2)
      IF(GRC.LT.1.E-10) GRC=1.E-10
      SC = 1.+CHISQC/CHISQ
      ST = ST/SC
      IF(CHISQC.GT.CHISQ)              GO TO 20
      COSI=(GRX*GRXC+GRY*GRYC)/GR/GRC
      ST = ST*SC/(1.4-COSI)
      X1 = XC
      Y1 = YC
      GRX=GRXC
      GRY=GRYC
      GR = GRC
      CHISQ=CHISQC
      IF(CHISQ.LT.CHSTOP)              GO TO 101
   20 STEP=ST*GR
      IF(STEP.LT.STPMIN)               GO TO 101
      IF(STEP.GT.1.)                   ST = 1./GR
      XC = X1 - ST*GRX
      YC = Y1 - ST*GRY
  100 CONTINUE
  101 CHISQ=CHISQ/DOF
      IF(CHISQ.GT.CSQCUT) CALL TWOGAM(NADC,IADC,CHISQ,E1,X1,Y1,E2,X2,Y2)
      RETURN
      END
C=======================================================================
      SUBROUTINE TWOGAM(NADC,IADC,CHISQ,E1,X1,Y1,E2,X2,Y2)
C   28 OCT.1981.     TWO GAMMA FIT PROGRAM.
CORRECTED 28 MAR.83.
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
      DIMENSION IADC(2,NADC)
      DATA STPMIN/.005/,CONST/1.5/,DELCH/6./,EMIN/20./,ST0/.00005/
      DATA CHMIN/1./
C
C     WRITE(5,7)
C   7 FORMAT(1X,'*** TWOGAM ***')
      IF(NADC.LE.3)                   RETURN
      DOF=NADC-5
      IF(DOF.LT.1.)                   DOF=1.
      CHSTOP=CHMIN*DOF
C   START POINT CHOOSING.
      E=0.
      XX=0.
      YY=0.
      DO 1 I=1,NADC
      A=IADC(1,I)
      IXY=IADC(2,I)
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      XX=XX + A*IX
      YY=YY + A*IY
      E = E + A
    1 CONTINUE
      XX=XY(XX/E)
      YY=XY(YY/E)
      EXX=0.
      EYY=0.
      EXY=0.
      DO 2 I=1,NADC
      A=IADC(1,I)
      IXY=IADC(2,I)
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      DX=IX-XX
      DY=IY-YY
      EXX=EXX + A*DX**2
      EYY=EYY + A*DY**2
      EXY=EXY + A*DX*DY
    2 CONTINUE
      SQR=SQRT(4.*EXY**2+(EXX-EYY)**2)
      EUU=(EXX+EYY+SQR)/2.
      COS2FI=1.
      IF(SQR.GT.1.E-10) COS2FI=(EXX-EYY)/SQR
      COSFI=SQRT((1.+COS2FI)/2.)
      SINFI=SQRT((1.-COS2FI)/2.)
      IF(EXY.LT.0.) SINFI=-SINFI
      EU3=0.
      DO 3 I=1,NADC
      A=IADC(1,I)
      IXY=IADC(2,I)
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      DX=IX-XX
      DY=IY-YY
      DU=DX*COSFI+DY*SINFI
    3 EU3=EU3+A*DU**3
      C=E*EU3**2/EUU**3/2.
      SIGN=1.
      IF(EU3.LT.0.) SIGN=-1.
      R=1.+C+SIGN*SQRT(2.*C+C**2)
      IF(ABS(R-1.).GT..1) U=EU3/EUU*(R+1.)/(R-1.)
      IF(ABS(R-1.).LE..1) U=SQRT(SQR/E/R)*(1.+R)
      E2C=E/(1.+R)
      E1C=E-E2C
      U1=-U/(1.+R)
      U2=U+U1
      X1C=XX+U1*COSFI
      Y1C=YY+U1*SINFI
      X2C=XX+U2*COSFI
      Y2C=YY+U2*SINFI
C
C   START OF ITERATIONS.
C
      ST = ST0
      CH = 1.E20
      GRX1 = 0.
      GRY1 = 0.
      GRX2 = 0.
      GRY2 = 0.
      GRE  = 0.
      GR   = 1.
      DO 100 ITER=1,50
      CHISQC= 0.
      GRX1C = 0.
      GRY1C = 0.
      GRX2C = 0.
      GRY2C = 0.
      GREC  = 0.
      DO 10 I=1,NADC
      IXY=IADC(2,I)
      IX=IXY/NcelY
      IY=IXY-IX*NcelY
      DX1 = X1C-IX
      DY1 = Y1C-IY
      CALL AG(E1C,DX1,DY1,A1,GX1,GY1)
      DX2 = X2C-IX
      DY2 = Y2C-IY
      CALL AG(E2C,DX2,DY2,A2,GX2,GY2)
      EX = IADC(1,I)
      A  = A1 + A2
      D  = CONST * A * (1.-A/E)
      IF(D.LT.0.) D = 0.
      D  =     9.
      DA = A - EX
      CHISQC=CHISQC+DA**2/D
      DD = DA/D
      DD = DD*(2.-DD*CONST*(1.-2.*A/E))
      GRX1C = GRX1C + GX1*DD
      GRY1C = GRY1C + GY1*DD
      GRX2C = GRX2C + GX2*DD
      GRY2C = GRY2C + GY2*DD
      GREC  = GREC  + (A1/E1C-A2/E2C)*E*DD
   10 CONTINUE
      GRC=SQRT(GRX1C**2+GRY1C**2+GRX2C**2+GRY2C**2+GREC**2)
      IF(GRC.LT.1.E-10) GRC=1.E-10
      SC = 1.+CHISQC/CH
      ST = ST/SC
      IF(CHISQC.GT.CH)              GO TO 20
      COSI=(GRX1*GRX1C+GRY1*GRY1C+GRX2*GRX2C+GRY2*GRY2C+GRE*GREC)/GR/GRC
      ST = ST*SC/(1.4-COSI)
      EE1 = E1C
      XX1 = X1C
      YY1 = Y1C
      EE2 = E2C
      XX2 = X2C
      YY2 = Y2C
      CH  = CHISQC
      IF(CH.LT.CHSTOP)               GO TO 101
      GRX1= GRX1C
      GRY1= GRY1C
      GRX2= GRX2C
      GRY2= GRY2C
      GRE = GREC
      GR  = GRC
   20 STEP= ST*GR
      IF(STEP.LT.STPMIN)            GO TO 101
      DE  = ST  * GRE * E
      E1C = EE1 - DE
      E2C = EE2 + DE
      IF(E1C.GT.EMIN.AND.E2C.GT.EMIN) GO TO 25
      ST  = ST / 2.
      GO TO 20
   25 X1C = XX1 - ST*GRX1
      Y1C = YY1 - ST*GRY1
      X2C = XX2 - ST*GRX2
      Y2C = YY2 - ST*GRY2
  100 CONTINUE
  101 IF(CH.GE.CHISQ*(NADC-2)-DELCH)         RETURN
      CHISQ=CH/DOF
      E1 = EE1
      X1 = XX1
      Y1 = YY1
      E2 = EE2
      X2 = XX2
      Y2 = YY2
      RETURN
      END
C=======================================================================
      FUNCTION XY(XI)
C   MADE FROM 25 GEV ELECTRON CALIBRATION WITHOUT HODOSCOPE.  RUN APR.1981.
      I=XI
      X=XI-I-.5
      XI2=X*X
      XY=X*(.22466+XI2*(3.74014+XI2*(-25.88467+
     + XI2*(166.90556-XI2*294.35077))))+I+.5
      RETURN
      END
C
C=======================================================================
      SUBROUTINE  AG(Ei,Xi,Yi,Ai,GXi,GYi)
      REAL*4         Ei,Xi,Yi,Ai,GXi,GYi  !...Interface...!
      REAL*8         E ,X ,Y ,A ,GX ,GY   !..Calculation..!
*  ..........................................................
*  ...It uses Lednev's Amplitude calculation function.Fcml..
*  ...To Calculate an Ammplitude (A) and Garadient (GX,GY)...
*  ...i use REAL*4 for inteface purpose only ................
*  ...calculations in REAL*8.................................
*  ...................typed by Sadovsky.A.S..................
      REAL*8 Fcml,Dx,Dy,D
      INTEGER i
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
C
      Dx= CelZ/2.
      Dy= CelY/2.
C
      X=Xi*CelZ  !....REAL*4....!
      Y=Yi*CelY  !.....to.......!
      E=Ei       !..............!
      A=Ai       !....REAL*8....!
      GX=GXi     !..Convertion..!
      GY=GYi     !..............!

      A = Fcml(x+dx,y+dy) - Fcml(x+dx,y-dy) -
     +    Fcml(x-dx,y+dy) + Fcml(x-dx,y-dy)
      Ai= A*E
C
      GX= GradX(x+dx,y+dy) - GradX(x+dx,y-dy) -
     +    GradX(x-dx,y+dy) + GradX(x-dx,y-dy)
      GXi=GX*E*E
C
      GY= GradY(x+dx,y+dy) - GradY(x+dx,y-dy) -
     +    GradY(x-dx,y+dy) + GradY(x-dx,y-dy)
      GYi=GY*E*E
C
      RETURN
      END
C=======================================================================
      FUNCTION Fcml(X,Y)
C                     *..Cumuliative function calculation in REAL*8..*
C
      REAL*8   Fcml,X,Y
      REAL*8    A,b,Fff
      DATA      A,b / 1.0, 1.082 /
C-                                !...Cumuliative function calculation.......!
C
C                                 !  Kumulyativnaya funkciya ne pravil'naya, !
C                                    no poluchaemoe pri etom enrgovydelenie
C                                          v yachejke PRAVIL'NOE !!!
C
      Fff  = A*DATAN(x*y/(b*DSQRT(b*b+x*x+y*y)))
      Fcml = Fff/6.2831853071796D0
      RETURN
      END
C=======================================================================
      FUNCTION GradX(X,Y)
C                      *..Cumuliative function Gradient_X in REAL*8..*
C
      REAL*8   GradX,X,Y
      REAL*8   a,b,skv,s,Gradient
      DATA     a,b / 1.0, 1.082 /
C
      skv      = b*b + x*x + y*y
      Gradient = a*y*(1.-x*x/skv)*b*DSQRT(skv)/(b*b*skv +x*x*y*y)
      GradX    = Gradient/6.2831853071796d0
      RETURN
      END
C=======================================================================
      FUNCTION GradY(X,Y)
C                      *..Cumuliative function Gradient_Y in REAL*8..*
C
C-    COMMON/CONSTS/ BEAM,ANOR,ANORT,DELTA,CELL,DSST,NCNX,NCNY,STCO
      REAL*8   GradY,X,Y
      REAL*8   a,b,skv,s,Gradient
      DATA     a,b / 1.0, 1.082 /
C
      skv      = b*b + x*x + y*y
      Gradient = a*x*(1.-y*y/skv)*b*DSQRT(skv)/(b*b*skv +x*x*y*y)
      GradY    = Gradient/6.2831853071796d0
      RETURN
      END
C=======================================================================
      FUNCTION AMP(Ei,Xi,Yi)
C!-
C     *..Amplitude calculation in REAL*4 using cumuliative function..*
C     *........to calculate amplitude only in gams cell..............*
C
      COMMON/CONSTA/ NcelZ,NcelY,IDELT,IDELA,CelZ,CelY
      REAL*4   AMP,Ei,Xi,Yi,Fcm
      REAL*8   Fcml,E,X ,Y, Dx,Dy,D
      INTEGER i
C
      Dx= CelZ/2.
      Dy= CelY/2.
      x = Xi*CelZ
      y = Yi*CelY
C
      E = Ei
      A = Fcml(x+dx,y+dy) - Fcml(x+dx,y-dy) -
     +    Fcml(x-dx,y+dy) + Fcml(x-dx,y-dy)
      AMP = A*E
      RETURN
      END

C=======================================================================
      SUBROUTINE VARDEL(NWGAMS)
      COMMON/DELTA/ IDLT(40000)
      COMMON/CLUSTER/ NCL,LENCL(5000)
      COMMON/EVENT/ IHDR(12),IBUF(30035)
      EQUIVALENCE(IHDR(3),NREC),(IHDR(5),NMT),(IHDR(6),IPHOD),(IHDR(7),
     ,IPSC),(IHDR(8),IPGS),(IHDR(9),IPGM),(IHDR(10),IPTG),(IBUF(3),NEV)
      DATA MINDLT/2/,NTRIG/1000/,LENMAX/6/
C
      IF(NWGAMS.LE.2.OR.NCL.LT.1)      RETURN
      ITRIG=ITRIG+1
      IF(ITRIG.LE.NTRIG)               GO TO 20
      ITRIG=0
      DO 10 I=1,1
      IF(IDLT(I).GT.MINDLT)            IDLT(I)=IDLT(I)-1
   10 CONTINUE
   20 IPCL=IPGM+2
      DO 50 ICL=1,NCL
      IF(LENCL(ICL).NE.1)              GO TO 30
      ICNT=IBUF(IPCL+1)+1
      IDLT(ICNT)=IDLT(ICNT)+1
      GO TO 50
   30 IF(LENCL(ICL).GT.LENMAX)         GO TO 50
      IPEND=IPCL+LENCL(ICL)*2-2
      DO 40 I=IPCL,IPEND,2
      ICNT=IBUF(I+1)+1
      IF(IDLT(ICNT).GT.MINDLT)         IDLT(ICNT)=IDLT(ICNT)-1
   40 CONTINUE
   50 IPCL=IPCL+LENCL(ICL)*2
      RETURN
      END
C=======================================================================

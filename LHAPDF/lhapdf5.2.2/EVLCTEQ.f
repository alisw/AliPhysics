      SUBROUTINE CtLhALFSET (QS, ALFS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL CtLhRTALF
      COMMON / CtLhRTALFC / ALFST, JORD, NEFF
      DATA ALAM, BLAM, ERR / 0.01, 10.0, 0.02 /
      QST   = QS
      ALFST = ALFS
      CALL CtLhParQcd (2, 'ORDR', ORDR, IR1)
      JORD  = ORDR
      NEFF = LhCtNFL(QS)
      EFLLN  = CtLhQZBRNT (CtLhRTALF, ALAM, BLAM, ERR, IR2)
      EFFLAM = QS / EXP (EFLLN)
      CALL CtLhSETL1 (NEFF, EFFLAM)
      END
      FUNCTION CtLhALPI (AMU)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
      DATA IW1, IW2 / 2*0 /
      IF(.NOT.SET) CALL CtLhLAMCWZ
      NEFF = LhCtNFL(AMU)
      ALM  = ALAM(NEFF)
      CtLhALPI = CtLhALPQCD (NORDER, NEFF, AMU/ALM, IRT)
      IF (IRT .EQ. 1) THEN
         CALL CtLhWARNR (IW1, 'AMU < ALAM in CtLhALPI', 'AMU', AMU,
     >              ALM, BIG, 1)
      ELSEIF (IRT .EQ. 2) THEN
         CALL CtLhWARNR (IW2, 'CtLhALPI > 3; Be aware!', 'CtLhALPI', 
     >  CtLhALPI, D0, D1, 0)
      ENDIF
      RETURN
      END
      FUNCTION CtLhALPQCD (IRDR, NF, RML, IRT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
      PARAMETER (CG = 3.0d0, TR = 0.5d0, CF = 4.0d0/3.0d0)
      IRT = 0
      IF (IRDR .LT. 1 .OR. IRDR .GT. 2) THEN
        print *,
     >  'Order out of range in CtLhALPQCD: IRDR = ', IRDR
        STOP
      ENDIF
      B0 = (11.d0*CG  - 2.* NF) / 3.d0
      B1 = (34.d0*CG**2 - 10.d0*CG*NF - 6.d0*CF*NF) / 3.d0
      RM2 = RML**2
      IF (RM2 .LE. 1.) THEN
         IRT = 1
         CtLhALPQCD = 99.
         RETURN
      ENDIF
      ALN = LOG (RM2)
      AL = 4.d0/ B0 / ALN
      IF (IRDR .GE. 2) AL = AL * (1.d0-B1*LOG(ALN) / ALN / B0**2)
      IF (AL .GE. 3.) THEN
         IRT = 2
      ENDIF
      CtLhALPQCD = AL
      RETURN
      END
      FUNCTION CtLhAMHATF(I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      IF (.NOT.SET) CALL CtLhLAMCWZ
      IF ((I.LE.0).OR.(I.GT.9)) THEN
         print *,'warning I OUT OF RANGE IN CtLhAMHATF'
         CtLhAMHATF = 0
      ELSE
         CtLhAMHATF = AMHAT(I)
      ENDIF
      RETURN
      END
      FUNCTION CtLhDXDZ (Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      DATA HUGE, IWRN / 1.E20, 0 /
      ZZ = Z
      X = CtLhXFRMZ (ZZ)
      TEM = CtLhDZDX (X)
      IF     (TEM .NE. D0) THEN
        TMP = D1 / TEM
      Else
      CALL CtLhWARNR(IWRN, 'CtLhDXDZ singular in CtLhDXDZ; set=HUGE',
     >             'Z', Z, D0, D1, 0)
        TMP = HUGE
      EndIf
      CtLhDXDZ = TMP
      RETURN
      END
      SUBROUTINE CtLhEVLPAR (IACT, NAME, VALUE, IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) NAME
      IRET = 1
      IF     (IACT .EQ. 0) THEN
              WRITE ( NINT(VALUE) , 101)
  101         FORMAT (/ ' Initiation parameters:   Qini, Ipd0, Ihdn ' /
     >                  ' Maximum Q, Order of Alpha:     Qmax, IKNL ' /
     >                  ' X- mesh parameters   :   Xmin, Xcr,   Nx  ' /
     >                  ' LnQ-mesh parameters  :         Nt,   Jt   ' /
     >                  ' # of parton flavors  :         NfMx       ' /)
              IRET = 4
      ElseIF (IACT .EQ. 1) THEN
              CALL CtLhEVLSET (NAME, VALUE, IRET)
      Else
	print *,'fatal evlpar'
	stop
      EndIf
      RETURN
      END
      SUBROUTINE CtLhEVLSET (NAME, VALUE, IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      CHARACTER*(*) NAME
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtQARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
     > / PdfSwh / Iset, IpdMod, Iptn0, NuIni
      IRET = 1
      IF     (NAME .EQ. 'QINI')  THEN
          IF (VALUE .LE. 0) GOTO 12
          QINI = VALUE
      ElseIF (NAME .EQ. 'IPD0')  THEN
          ITEM = NINT(VALUE)
          IF (Item .Eq. 10 .or. Item .Eq. 11) GOTO 12
          IPD0 = ITEM
      ElseIF (NAME .EQ. 'IHDN') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. -1 .OR. ITEM .GT. 5) GOTO 12
          IHDN = ITEM
      ElseIF (NAME .EQ. 'QMAX')  THEN
          IF (VALUE .LE. QINI) GOTO 12
          QMAX = VALUE
      ElseIF (NAME .EQ. 'IKNL') THEN
          ITMP = NINT(VALUE)
          ITEM = ABS(ITMP)
          IF (ITEM.NE.1.AND.ITEM.NE.2) GOTO 12
          IKNL = ITMP
      ElseIF (NAME .EQ. 'XCR') THEN
          IF (VALUE .LT. XMIN .OR. VALUE .GT. 10.) GOTO 12
          XCR = VALUE
          LSTX = .FALSE.
      ElseIF (NAME .EQ. 'XMIN') THEN
          IF (VALUE .LT. 1D-7 .OR. VALUE .GT. 1D0) GOTO 12
          XMIN = VALUE
          LSTX = .FALSE.
      ElseIF (NAME .EQ. 'NX') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 10 .OR. ITEM .GT. MXX-1) GOTO 12
          NX = ITEM
          LSTX = .FALSE.
      ElseIF (NAME .EQ. 'NT') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 2 .OR. ITEM .GT. MXQ) GOTO 12
          NT = ITEM
      ElseIF (NAME .EQ. 'JT') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 1 .OR. ITEM .GT. 5) GOTO 12
          JT = ITEM
      ElseIF (NAME .EQ. 'NFMX') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 1 .OR. ITEM .GT. MXPN) GOTO 12
          NfMx = ITEM
      ElseIF (NAME .EQ. 'IPDMOD') THEN
          ITEM = NINT(VALUE)
          IF (Abs(Item) .Gt. 1) GOTO 12
          IpdMod = ITEM
      ElseIF (NAME .EQ. 'IPTN0') THEN
          ITEM = NINT(VALUE)
          IF (ABS(ITEM) .GT. MXF) GOTO 12
          IPTN0 = ITEM
      ElseIF (NAME .EQ. 'NUINI') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LE. 0) GOTO 12
          NuIni = ITEM
      Else
          IRET = 0
      EndIf
      RETURN
   12 IRET = 2
      RETURN
      END
      SUBROUTINE CtLhEVOLVE (FINI, IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      include 'parmsetup.inc'
      LOGICAL LSTX
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtQARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / LhCtQARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / LhCtPEVLDT / UPD(MXPQX,nmxset), KF, Nelmt
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtVARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
      DIMENSION QRKP(MXF)
      DIMENSION JI(-MXF : MXF+1)
      EXTERNAL LhCtNSRHSP, LhCtNSRHSM, FINI
      DATA DZER / 0.0 /
	save nxsave, ntsave, jtsave, ngsave, 
     &	     xcrsave, xminsave, qinisave, qmaxsave, ientry, ishow
	data ientry / 0 /
	dimension xvsave(0:MXX), qvsave(0:MXQ), tvsave(0:MXQ)
c
       call getnset(iset)
c	
	ientry = ientry + 1
	if(ientry .lt. 100) then
666	   format(1x,'enter evolve',i3)
	elseif(ientry .eq. 100) then
667	   format(1x,'enter evolve',i3,' further suppressed')
	endif
	ishow = 0		!set to no display
	if(ientry .eq. 1) then
	   ishow = 1		!turn display on
	endif
   11 IRET = 0
      IF (IHDN .LE. 4) THEN
        MXVAL = 2
      ElseIF (IHDN .LE. 6) THEN
        MXVAL = 3
      EndIf
      IF (.NOT. LSTX) CALL CtLhXARRAY
      DLX = 1.D0 / NX
      J10 = NX / 10
      CALL CtLhPARPDF (2, 'ALAM', AL, IR)
      CALL CtLhQARRAY (NINI)
      NFSN = NFMX + 1
      KF = 2 * NFMX + 2
      Nelmt = KF * (Nt+1) * (Nx+1)
      DO 101 IFLV = -NFMX, NFMX+1
        JFL = NFMX + IFLV
        JI(IFLV) = JFL * (NT+1) * (NX+1)
  101 CONTINUE
    3 DO 31 IZ = 1, NX
        UPD(JI(0)+IZ+1,iset) = FINI (0, XV(IZ))
        UPD(JI(NFSN)+IZ+1,iset) = 0
        IF (NFMX .EQ. 0) GOTO 31
        DO 331 IFLV = 1, NINI
          A = FINI ( IFLV, XV(IZ))
          B = FINI (-IFLV, XV(IZ))
          QRKP (IFLV) = A + B
          UPD(JI(NFSN)+IZ+1,iset) = 
     >       UPD(JI(NFSN)+IZ+1,iset) + QRKP (IFLV)
          UPD(JI(-IFLV)+IZ+1,iset) = A - B
  331   CONTINUE
        DO 332 IFLV = 1, NINI
           UPD(JI( IFLV)+IZ+1,iset) = 
     >        QRKP(IFLV) - UPD(JI(NFSN)+IZ+1,iset)/NINI
  332   CONTINUE
   31 CONTINUE
      DO 21 NEFF = NINI, NFMX
          IF (IKNL .EQ. 2) CALL CtLhSTUPKL (NEFF)
          ICNT = NEFF - NINI + 1
          IF (NTN(ICNT) .EQ. 0) GOTO 21
          NITR = NTN (ICNT)
          DT   = DTN (ICNT)
          TIN  = TLN (ICNT)
          CALL CtLhSNEVL (IKNL, NX, NITR, JT, DT, TIN, NEFF,
     >    UPD(JI(NFSN)+2,iset), UPD(JI(0)+2,iset),
     >    UPD(JI(NFSN)+1,iset), UPD(JI(0)+1,iset))
          IF (NEFF .EQ. 0) GOTO 88
    5     DO 333 IFLV = 1, NEFF
           CALL CtLhNSEVL (LhCtNSRHSP, IKNL, NX, NITR, JT, DT, TIN, 
     >     NEFF, UPD(JI( IFLV)+2,iset), UPD(JI( IFLV)+1,iset))
           IF (IFLV .LE. MXVAL)
     >     CALL CtLhNSEVL (LhCtNSRHSM, IKNL, NX, NITR, JT, DT, TIN, 
     >     NEFF, UPD(JI(-IFLV)+2,iset), UPD(JI(-IFLV)+1,iset))
           DO 55 IS = 0, NITR
           DO 56 IX = 0, NX
             TP = UPD (IS*(NX+1) + IX + 1 + JI( IFLV),iset)
             TS = UPD (IS*(NX+1) + IX + 1 + JI( NFSN),iset) / NEFF
             TP = TP + TS
             IF (IKNL .GT. 0) TP = MAX (TP, DZER)
             IF (IFLV .LE. MXVAL) THEN
                TM = UPD (IS*(NX+1) + IX + 1 + JI(-IFLV),iset)
                IF (IKNL .GT. 0) THEN
                  TM = MAX (TM, DZER)
                  TP = MAX (TP, TM)
                EndIf
             Else
                TM = 0.
             EndIf
             UPD (JI( IFLV) + IS*(NX+1) + IX + 1,iset) = (TP + TM)/2.
             UPD (JI(-IFLV) + IS*(NX+1) + IX + 1,iset) = (TP - TM)/2.
   56      CONTINUE
   55      CONTINUE
333      CONTINUE
        DO 334 IFLV = NEFF + 1, NFMX
          DO 57 IS = 0, NITR
          DO 58 IX = 0, NX
            UPD(JI( IFLV) + IS*(NX+1) + IX + 1,iset) = 0
            UPD(JI(-IFLV) + IS*(NX+1) + IX + 1,iset) = 0
   58     CONTINUE
   57     CONTINUE
  334   CONTINUE
   88   CONTINUE
        IF (NFMX .EQ. NEFF) GOTO 21
        DO 335 IFLV = -NFMX, NFMX+1
           JI(IFLV) = JI(IFLV) + NITR * (NX+1)
  335   CONTINUE
        CALL CtLhHQRK (NX, TT, NEFF+1, UPD(JI(0)+2,iset),
     >     UPD(JI(NEFF+1)+2,iset))
        DO 32 IZ = 1, NX
         QRKP (NEFF+1) = 2. * UPD(JI( NEFF+1) + IZ + 1,iset)
         UPD (JI(NFSN)+IZ+1,iset) = UPD (JI(NFSN)+IZ+1,iset)
     >        + QRKP (NEFF+1)
         VS00 =  UPD (JI(NFSN)+IZ+1,iset) / (NEFF+1)
         UPD(JI( NEFF+1) + IZ + 1,iset) = QRKP(NEFF+1) - VS00
         DO 321 IFL = 1, NEFF
           A = UPD(JI( IFL)+IZ+1,iset)
           B = UPD(JI(-IFL)+IZ+1,iset)
           QRKP(IFL) = A + B
           UPD(JI( IFL)+IZ+1,iset) = QRKP(IFL) - VS00
           IF (IFL .LE. MXVAL)  UPD(JI(-IFL)+IZ+1,iset) = A - B
  321    CONTINUE
   32   CONTINUE
   21 CONTINUE
	if(ientry .eq. 1) then
	   nxsave = nx
	   ntsave = nt
	   jtsave = jt
	   ngsave = ng
	   xcrsave = xcr
	   xminsave = xmin
	   qinisave = qini
	   qmaxsave = qmax
	endif
	if((nx .ne. nxsave) .or.
     &	   (nt .ne. ntsave) .or.
     &	   (jt .ne. jtsave) .or.
     &	   (ng .ne. ngsave) .or.
     &	   (xcr .ne. xcrsave) .or.
     &	   (xmin .ne. xminsave) .or.
     &	   (qini .ne. qinisave) .or.
     &	   (qmax .ne. qmaxsave)) then
	   write(6,669) nx, nt, jt, ng, xcr, xmin, 
     &	                qini, qmax, ientry
669	   format(1x,'evolve.f:  nx,nt,jt,ng=',4i4,
     &	             ' xcr,xmin=',2f9.6,
     &	             ' qini, qmax',f7.4,1x,e12.5,' ientry=',i6)
	   nxsave = nx
	   ntsave = nt
	   jtsave = jt
	   ngsave = ng
	   qinisave = qini
	   qmaxsave = qmax
	   xcrsave = xcr
	   xminsave = xmin
	endif
	ixshow = 0
	do i = 0, mxx
	   if(xvsave(i) .ne. xv(i)) then
	      ixshow = 1
	      xvsave(i) = xv(i)
	   endif
	enddo
	iqshow = 0
	itshow = 0
	do i = 0, mxq
	   if(qvsave(i) .ne. qv(i)) then
	      iqshow = 1
	      qvsave(i) = qv(i)
	   endif
	   if(tvsave(i) .ne. tv(i)) then
	      itshow = 1
	      tvsave(i) = tv(i)
	   endif
	enddo
      Return
      End
      FUNCTION CtLhFINTRP (FF,  X0, DX, NX,  XV,  ERR, IR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 3)
      DIMENSION FF (0:NX), XX(MX)
      DATA SML, XX / 1.D-5,  0., 1.0, 2.0 /
      DATA  IW1, IW3, IW5 / 3 * 0 /
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      CtLhFINTRP = 0.
      IF (NX .LT. 1) THEN
         CALL CtLhWARNI(IW1, 'Nx < 1, error in CtLhFINTRP.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN(NX+1, MX)
      ENDIF
      IF (DX .LE. 0) THEN
         CALL CtLhWARNR(IW3, 'DX < 0, error in CtLhFINTRP.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
      XM = X0 + DX * NX
      IF (X .LT. X0-SML .OR. X .GT. XM+SML) THEN
        CALL CtLhWARNR(IW5,
     >     'X out of range in CtLhFINTRP, Extrapolation used.',
     >     'X',X,X0,XM,1)
      IR = 3
      ENDIF
      TX = (X - X0) / DX
      IF (TX .LE. 1.) THEN
        IX = 0
      ELSEIF (TX .GE. ANX-1.) THEN
        IX = NX - 2
      ELSE
        IX = TX
      ENDIF
      DDX = TX - IX
      CALL CtLhRATINT (XX, FF(IX), MNX, DDX, TEM, ERR)
      CtLhFINTRP = TEM
      RETURN
      END
      FUNCTION CtLhGausInt(F,XL,XR,AERR,RERR,ERR,IRT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        DIMENSION XLIMS(100), R(93), W(93)
        INTEGER PTR(4),NORD(4)
        external f
        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860,.6612093865,.9324695142,
     1 .1252334085,.3678314990,.5873179543,.7699026742,.9041172563,
     1 .9815606342,.0640568929,.1911188675,.3150426797,.4337935076,
     1 .5454214714,.6480936519,.7401241916,.8200019860,.8864155270,
     1 .9382745520,.9747285560,.9951872200,.0323801710,.0970046992,
     1 .1612223561,.2247637903,.2873624873,.3487558863,.4086864820,
     1 .4669029048,.5231609747,.5772247261,.6288673968,.6778723796,
     1 .7240341309,.7671590325,.8070662040,.8435882616,.8765720203,
     1 .9058791367,.9313866907,.9529877032,.9705915925,.9841245837,
     1 .9935301723,.9987710073,.0162767488,.0488129851,.0812974955,
     1 .1136958501,.1459737146,.1780968824,.2100313105,.2417431561,
     1 .2731988126,.3043649444,.3352085229,.3656968614,.3957976498,
     1 .4254789884,.4547094222,.4834579739,.5116941772,.5393881083,
     1 .5665104186,.5930323648,.6189258401,.6441634037,.6687183100,
     1 .6925645366,.7156768123,.7380306437,.7596023411,.7803690438,
     1 .8003087441,.8194003107,.8376235112,.8549590334,.8713885059,
     1 .8868945174,.9014606353,.9150714231,.9277124567,.9393703398,
     1 .9500327178,.9596882914,.9683268285,.9759391746,.9825172636,
     1 .9880541263,.9925439003,.9959818430,.9983643759,.9996895039/
        DATA W/.4679139346,.3607615730,.1713244924,
     1 .2491470458,.2334925365,.2031674267,.1600783285,.1069393260,
     1 .0471753364,.1279381953,.1258374563,.1216704729,.1155056681,
     1 .1074442701,.0976186521,.0861901615,.0733464814,.0592985849,
     1 .0442774388,.0285313886,.0123412298,.0647376968,.0644661644,
     1 .0639242386,.0631141923,.0620394232,.0607044392,.0591148397,
     1 .0572772921,.0551995037,.0528901894,.0503590356,.0476166585,
     1 .0446745609,.0415450829,.0382413511,.0347772226,.0311672278,
     1 .0274265097,.0235707608,.0196161605,.0155793157,.0114772346,
     1 .0073275539,.0031533461,.0325506145,.0325161187,.0324471637,
     1 .0323438226,.0322062048,.0320344562,.0318287589,.0315893308,
     1 .0313164256,.0310103326,.0306713761,.0302999154,.0298963441,
     1 .0294610900,.0289946142,.0284974111,.0279700076,.0274129627,
     1 .0268268667,.0262123407,.0255700360,.0249006332,.0242048418,
     1 .0234833991,.0227370697,.0219666444,.0211729399,.0203567972,
     1 .0195190811,.0186606796,.0177825023,.0168854799,.0159705629,
     1 .0150387210,.0140909418,.0131282296,.0121516047,.0111621020,
     1 .0101607705,.0091486712,.0081268769,.0070964708,.0060585455,
     1 .0050142027,.0039645543,.0029107318,.0018539608,.0007967921/
        DATA TOLABS,TOLREL,NMAX/1.E-35,5.E-4,100/
        TOLABS=AERR
        TOLREL=RERR
     
        CtLhGausInt=0.
        NLIMS=2
        XLIMS(1)=XL
        XLIMS(2)=XR
10      AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2D0
        BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2D0
        TVAL=0.
        DO 15 I=1,3
15      TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        TVAL=TVAL*AA
        DO 25 J=1,4
        VAL=0.
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
20      VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
                CtLhGausInt=CtLhGausInt+VAL
                NLIMS=NLIMS-2
                IF (NLIMS.NE.0) GO TO 10
                RETURN
                END IF
25      TVAL=VAL
        IF (NMAX.EQ.2) THEN
                CtLhGausInt=VAL
                RETURN
                END IF
        IF (NLIMS.GT.(NMAX-2)) THEN
                write(*,50) CtLhGausInt,NMAX,BB-AA,BB+AA
                RETURN
                END IF
        XLIMS(NLIMS+1)=BB
        XLIMS(NLIMS+2)=BB+AA
        XLIMS(NLIMS)=BB
        NLIMS=NLIMS+2
        GO TO 10
50      FORMAT (' CtLhGausInt FAILS, CtLhGausInt,NMAX,XL,XR=',
     >            G15.7,I5,2G15.7)
        END
      SUBROUTINE CtLhHINTEG (NX, F, H)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtHINTEC / GH(NDG, MXX)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      DIMENSION F(NX), H(NX), G(MXX)
      DZ = 1D0 / (NX-1)
      DO 20 I = 1, NX-2
         NP = NX - I + 1
         TEM = GH(1,I)*F(I) + GH(2,I)*F(I+1) + GH(3,I)*F(I+2)
         DO 30 KZ = 3, NP
            IY = I + KZ - 1
            W = XA(I,1) / XA(IY,1)
            G(KZ) = DXTZ(IY)*(F(IY)-W*F(I))/(1.-W)
   30    CONTINUE
         HTEM = CtLhSMPSNA (NP-2, DZ, G(3), ERR)
         TEM1 = F(I) * ELY(I)
         H(I) = TEM + HTEM + TEM1
   20 CONTINUE
      H(NX-1) = F(NX) - F(NX-1) + F(NX-1) * (ELY(NX-1) - XA(NX-1,0))
      H(NX)   = 0
      RETURN
      END
      SUBROUTINE CtLhHQRK (NX, TT, NQRK, Y, F)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      DIMENSION Y(NX), F(NX)
      IF (NX .GT. 1) GOTO 11
   11 CONTINUE
      DO 230 IZ = 1, NX
        IF (NX .GT. 1) THEN
        F(IZ) = 0
        GOTO 230
        EndIf
  230 CONTINUE
      RETURN
      END
      SUBROUTINE CtLhINTEGR (NX, M, F,   G, IR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER MSG*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtVARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
      DIMENSION   F(NX), G(NX)
      DATA IWRN1, IWRN2 / 0, 0 /
      IRR = 0
      IF (NX .LT. 1 .OR. XA(NX-1,1) .EQ. 0D0) THEN
        MSG = 'NX out of range in CtLhINTEGR call'
        CALL CtLhWARNI (IWRN1, MSG, 'NX', NX, 0, MXX, 0)
        IRR = 1
      EndIf
      IF (M .LT. M1 .OR. M .GT. M2) THEN
        MSG ='Exponent M out of range in CtLhINTEGR'
        CALL CtLhWARNI (IWRN2, MSG, 'M', M, M1, M2, 1)
        IRR = 2
      EndIf
      G(NX) = 0D0
      TEM = H(1, NX-1, -M) * F(NX-2) + H(2, NX-1, -M) * F(NX-1)
     >    + H(3, NX-1, -M) * F(NX)
      IF (M .EQ. 0) THEN
         G(NX-1) = TEM
      Else
         G(NX-1) = TEM * XA(NX-1, M)
      EndIf
      DO 10 I = NX-2, 2, -1
         TEM = TEM + H(1,I,-M)*F(I-1) + H(2,I,-M)*F(I)
     >             + H(3,I,-M)*F(I+1) + H(4,I,-M)*F(I+2)
         IF (M .EQ. 0) THEN
            G(I) = TEM
         Else
            G(I) = TEM * XA(I, M)
         EndIf
   10 CONTINUE
      TEM = TEM + H(2,1,-M)*F(1) + H(3,1,-M)*F(2) + H(4,1,-M)*F(3)
      IF (M .EQ. 0) THEN
         G(1) = TEM
      Else
         G(1) = TEM * XA(1, M)
      EndIf
      IR = IRR
      RETURN
      END
      SUBROUTINE CtLhKERNEL
     >(XX, FF1, FG1, GF1, GG1, PNSP, PNSM, FF2, FG2, GF2, GG2, NFL, IRT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589793d0, PI2 = PI ** 2)
      PARAMETER (D0 = 0.0, D1 = 1.0)
      DATA CF, CG, TR, IWRN / 1.33333333333333d0, 3.0d0, 0.5d0, 0 /
      IRT = 0
      TRNF = TR * NFL
      X = XX
      IF (X .LE. 0. .OR. X .GE. 1.) THEN
        CALL CtLhWARNR(IWRN, 'X out of range in CtLhKERNEL', 'X', X,
     >             D0, D1, 1)
        IRT = 1
        RETURN
      EndIf
      XI = 1./ X
      X2 = X ** 2
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      SPEN2 = CtLhSPENC2 (X)
      FFP = (1.+ X2) * XM1I
      FGP = (2.- 2.* X + X2) / X
      GFP = 1. - 2.* X + 2.* X2
      GGP = XM1I + XI - 2. + X - X2
      FFM = (1.+ X2) * XP1I
      FGM = - (2.+ 2.* X + X2) / X
      GFM = 1. + 2.* X + 2.* X2
      GGM = XP1I - XI - 2. - X - X2
      FF1 = CF * FFP * (1.- X)
      FG1 = CF * FGP * X
      GF1 = 2.* TRNF * GFP
      GG1 = 2.* CG * GGP * X * (1.-X)
      PCF2 = -2.* FFP *XLN*XLN1M - (3.*XM1I + 2.*X)*XLN
     >     - (1.+X)/2.*XLN2 - 5.*(1.-X)
      PCFG = FFP * (XLN2 + 11.*XLN/3.+ 67./9.- PI**2 / 3.)
     >     + 2.*(1.+X) * XLN + 40.* (1.-X) / 3.
      PCFT = (FFP * (- XLN - 5./3.) - 2.*(1.-X)) * 2./ 3.
      PQQB = 2.* FFM * SPEN2 + 2.*(1.+X)*XLN + 4.*(1.-X)
      PQQB = (CF**2-CF*CG/2.) * PQQB
      PQQ2 = CF**2 * PCF2 + CF*CG * PCFG / 2. + CF*TRNF * PCFT
      PNSP = (PQQ2 + PQQB) * (1.-X)
      PNSM = (PQQ2 - PQQB) * (1.-X)
      FFCF2 = - 1. + X + (1.- 3.*X) * XLN / 2. - (1.+ X) * XLN2 / 2.
     >      - FFP * (3.* XLN / 2. + 2.* XLN * XLN1M)
     >      + FFM * 2.* SPEN2
      FFCFG = 14./3.* (1.-X)
     >      + FFP * (11./6.* XLN + XLN2 / 2. + 67./18. - PI2 / 6.)
     >      - FFM * SPEN2
      FFCFT = - 16./3. + 40./3.* X + (10.* X + 16./3.* X2 + 2.) * XLN
     >                 - 112./9.* X2 + 40./9./X - 2.* (1.+ X) * XLN2
     >      - FFP * (10./9. + 2./3. * XLN)
      FGCF2 = - 5./2.- 7./2.* X + (2.+ 7./2.* X) * XLN + (X/2.-1.)*XLN2
     >               - 2.* X * XLN1M
     >      - FGP * (3.* XLN1M + XLN1M ** 2)
      FGCFG = 28./9. + 65./18.* X + 44./9. * X2 - (12.+ 5.*X + 8./3.*X2)
     >                      * XLN + (4.+ X) * XLN2 + 2.* X * XLN1M
     >      + FGP * (-2.*XLN*XLN1M + XLN2/2. + 11./3.*XLN1M + XLN1M**2
     >               - PI2/6. + 0.5)
     >      + FGM * SPEN2
      FGCFT = -4./3.* X - FGP * (20./9.+ 4./3.*XLN1M)
      GFCFT = 4.- 9.*X + (-1.+ 4.*X)*XLN + (-1.+ 2.*X)*XLN2 + 4.*XLN1M
     >      + GFP * (-4.*XLN*XLN1M + 4.*XLN + 2.*XLN2 - 4.*XLN1M
     >               + 2.*XLN1M**2 - 2./3.* PI2 + 10.)
      GFCGT = 182./9.+ 14./9.*X + 40./9./X + (136./3.*X - 38./3.)*XLN
     >               - 4.*XLN1M - (2.+ 8.*X)*XLN2
     >      + GFP * (-XLN2 + 44./3.*XLN - 2.*XLN1M**2 + 4.*XLN1M
     >               + PI2/3. - 218./9.)
     >      + GFM * 2. * SPEN2
      GGCFT = -16.+ 8.*X + 20./3.*X2 + 4./3./X + (-6.-10.*X)*XLN
     >        - 2.* (1.+ X) * XLN2
      GGCGT = 2.- 2.*X + 26./9.*X2 - 26./9./X - 4./3.*(1.+X)*XLN
     >      - GGP * 20./9.
      GGCG2 = 27./2.*(1.-X) + 67./9.*(X2-XI) + 4.*(1.+X)*XLN2
     >              + (-25.+ 11.*X - 44.*X2)/3.*XLN
     >      + GGP * (67./9.- 4.*XLN*XLN1M + XLN2 - PI2/3.)
     >      + GGM * 2.* SPEN2
      FF2 = CF * TRNF * FFCFT + CF ** 2 * FFCF2 + CF * CG   * FFCFG
      FG2 = CF * TRNF * FGCFT + CF ** 2 * FGCF2 + CF * CG   * FGCFG
      GF2 = CF * TRNF * GFCFT                   + CG * TRNF * GFCGT
      GG2 = CF * TRNF * GGCFT + CG ** 2 * GGCG2 + CG * TRNF * GGCGT
      XLG = (LOG(1./(1.-X)) + 1.)
      XG2 = XLG ** 2
      FF2 = FF2 * X * (1.- X)
      FG2 = FG2 * X / XG2
      GF2 = GF2 * X / XG2
      GG2 = GG2 * X * (1.- X)
      RETURN
      END
      SUBROUTINE CtLhLAMCWZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      CALL CtLhSETL1 (NF, AL)
      END
      FUNCTION LhCtNAMQCD(NNAME)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER NNAME*(*), NAME*8
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
      NAME = NNAME
      LhCtNAMQCD=0
      IF ( (NAME .EQ. 'ALAM') .OR. (NAME .EQ. 'LAMB') .OR.
     1        (NAME .EQ. 'LAM') .OR. (NAME .EQ. 'LAMBDA') )
     2             LhCtNAMQCD=1
      IF ( (NAME .EQ. 'NFL') .OR. (NAME(1:3) .EQ. '#FL') .OR.
     1        (NAME .EQ. '# FL') )
     2             LhCtNAMQCD=2
      DO 10 I=1, 9
         IF (NAME .EQ. 'M'//CHAR(I+IASC0))
     1             LhCtNAMQCD=I+2
10       CONTINUE
      DO 20 I= 0, NF
         IF (NAME .EQ. 'LAM'//CHAR(I+IASC0))
     1             LhCtNAMQCD=I+13
20       CONTINUE
      IF (NAME(:3).EQ.'ORD' .OR. NAME(:3).EQ.'NRD') LhCtNAMQCD = 24
      RETURN
      END
      FUNCTION LhCtNFL(AMU)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      IF (.NOT. SET) CALL CtLhLAMCWZ
      LhCtNFL = NF - NHQ
      IF ((LhCtNFL .EQ. NF) .OR. (AMU .LE. AMN)) GOTO 20
      DO 10 I = NF - NHQ + 1, NF
         IF (AMU .GE. AMHAT(I)) THEN
            LhCtNFL = I
         ELSE
            GOTO 20
         ENDIF
10       CONTINUE
20    RETURN
      END
      SUBROUTINE CtLhNSEVL (RHS, IKNL,NX,NT,JT,DT,TIN,NEFF,U0,UN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      DIMENSION  U0(NX), UN(0:NX, 0:NT)
      DIMENSION  Y0(MXX), Y1(MXX), YP(MXX), F0(MXX), F1(MXX), FP(MXX)
      external rhs
      DDT = DT / JT
      IF (NX .GT. MXX) THEN
      WRITE (*,*) 'Nx =', NX, ' greater than Max pts in CtLhNSEVL.'
      STOP 'Program stopped in CtLhNSEVL'
      EndIf
      TMD = TIN + DT * NT / 2.
      AMU = EXP(TMD)
      TEM = 6./ (33.- 2.* NEFF) / CtLhALPI(AMU)
      TLAM = TMD - TEM
      DO 9 IX = 1, NX
      UN(IX, 0)  = U0(IX)
    9 CONTINUE
      UN(0, 0) = 3D0*U0(1) - 3D0*U0(2) - U0(1)
      TT = TIN
      DO 10 IZ = 1, NX
      Y0(IZ)   = U0(IZ)
   10 CONTINUE
      DO 20 IS = 1, NT
         DO 202 JS = 1, JT
            IRND = (IS-1) * JT + JS
            IF (IRND .EQ. 1) THEN
                CALL RHS (TT, Neff, Y0, F0)
                DO 250 IZ = 1, NX
                   Y0(IZ) = Y0(IZ) + DDT * F0(IZ)
  250           CONTINUE
                TT = TT + DDT
                CALL RHS (TT, NEFF, Y0, F1)
                DO 251 IZ = 1, NX
                   Y1(IZ) = U0(IZ) + DDT * (F0(IZ) + F1(IZ)) / 2D0
  251           CONTINUE
            Else
                CALL RHS (TT, NEFF, Y1, F1)
                DO 252 IZ = 1, NX
                   YP(IZ) = Y1(IZ) + DDT * (3D0 * F1(IZ) - F0(IZ)) / 2D0
  252           CONTINUE
                TT = TT + DDT
                CALL RHS (TT, NEFF, YP, FP)
                DO 253 IZ = 1, NX
                   Y1(IZ) = Y1(IZ) + DDT * (FP(IZ) + F1(IZ)) / 2D0
                   F0(IZ) = F1(IZ)
  253           CONTINUE
            EndIf
  202    CONTINUE
         DO 260 IZ = 1, NX
            UN (IZ, IS) = Y1(IZ)
  260    CONTINUE
         UN(0, IS) = 3D0*Y1(1) - 3D0*Y1(2) + Y1(3)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE LhCtNSRHSM (TT, NEFF, FI, FO)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtXYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / LhCtKRNL01 / AFF2(MXX),AFG2(MXX),AGF2(MXX),AGG2(MXX),
     >                  ANSP (MXX), ANSM (MXX), ZFG2, ZGF2, ZQQB
      COMMON / LhCtKRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      DIMENSION G1(MXX), FI(NX), FO(NX)
      DIMENSION W0(MXX), W1(MXX), WH(MXX)
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = CtLhALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
      CALL CtLhINTEGR (NX, 0, FI, W0, IR1)
      CALL CtLhINTEGR (NX, 1, FI, W1, IR2)
      CALL CtLhHINTEG (NX,    FI, WH)
      DO 230 IZ = 1, NX
      FO(IZ) = 2.* FI(IZ) + 4./3.* ( 2.* WH(IZ) - W0(IZ) - W1(IZ))
      FO(IZ) = CPL * FO(IZ)
  230 CONTINUE
      IF (IKNL .EQ. 2) THEN
      DZ = 1./ (NX - 1)
      DO 21 IX = 1, NX-1
        X = XV(IX)
        NP = NX - IX + 1
        IS = NP
        DO 31 KZ = 2, NP
          IY = IX + KZ - 1
          IT = NX - IY + 1
          XY = ZZ (IS, IT)
          G1(KZ) = PNS (IS,IT) * (FI(IY) - XY * FI(IX))
   31   CONTINUE
        TEM1 = CtLhSMPNOL (NP, DZ, G1, ERR)
        TMP2 = (TEM1 - FI(IX) * ANSM(IX)) * CPL2
        FO(IX) = FO(IX) + TMP2
   21 CONTINUE
      EndIf
      RETURN
      END
      SUBROUTINE LhCtNSRHSP (TT, NEFF, FI, FO)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtXYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / LhCtKRNL01 / AFF2(MXX),AFG2(MXX),AGF2(MXX),AGG2(MXX),
     >                  ANSP (MXX), ANSM (MXX), ZFG2, ZGF2, ZQQB
      COMMON / LhCtKRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      DIMENSION G1(MXX), FI(NX), FO(NX)
      DIMENSION W0(MXX), W1(MXX), WH(MXX)
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = CtLhALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
      CALL CtLhINTEGR (NX, 0, FI, W0, IR1)
      CALL CtLhINTEGR (NX, 1, FI, W1, IR2)
      CALL CtLhHINTEG (NX,    FI, WH)
      DO 230 IZ = 1, NX
      FO(IZ) = 2.* FI(IZ) + 4./3.* ( 2.* WH(IZ) - W0(IZ) - W1(IZ))
      FO(IZ) = CPL * FO(IZ)
  230 CONTINUE
      IF (IKNL .EQ. 2) THEN
      DZ = 1./ (NX - 1)
      DO 21 IX = 1, NX-1
        X = XV(IX)
        NP = NX - IX + 1
        DO 31 KZ = 2, NP
          IY = IX + KZ - 1
          XY = ZZ (NX-IX+1, NX-IY+1)
          G1(KZ) = PNS (IX,IY) * (FI(IY) - XY * FI(IX))
   31   CONTINUE
        TEM1 = CtLhSMPNOL (NP, DZ, G1, ERR)
        TMP2 = (TEM1 + FI(IX) * (-ANSP(IX) + ZQQB)) * CPL2
        FO(IX) = FO(IX) + TMP2
   21 CONTINUE
      EndIf
      RETURN
      END
      FUNCTION CtLhPARDIS (IPRTN, XX, QQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      include 'parmsetup.inc'
      Character Msg*80
      LOGICAL LSTX
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      PARAMETER (Smll = 1D-9)
	parameter(nqvec = 4)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtXYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtQARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / LhCtQARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / LhCtPEVLDT / UPD(MXPQX,nmxset), KF, Nelmt
      COMMON / LhCtCOMQMS / VALQMS(9)
      dimension fvec(4), fij(4)
      dimension xvpow(0:mxx)
      Data Iwrn1, Iwrn2, Iwrn3, OneP / 3*0, 1.00001 /
      data xpow / 0.3d0 /	!**** choice of interpolation variable
      data nxsave / 0 /
	save xvpow, nxsave
 	save xlast, qlast
 	save jq, jx, JLx, JLq, ss, sy2, sy3, s23, ty2, ty3
 	save const1 , const2, const3, const4, const5, const6
 	save tt, t13, t12, t23, t34 , t24, tmp1, tmp2, tdet
c
      call getnset(iset)
c
      if(nx .ne. nxsave) then
         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
	nxsave = nx
      endif

      X = XX
      Q = QQ

c enforce threshold early to improve speed...
	ii = iabs(IPRTN)
	if(ii .ne. 0) then
	   if(QQ .lt. VALQMS(ii) ) then
	      ctlhpardis = 0.d0
	      return
	   endif
	endif

c force pardis = 0.0d0 at exactly =1.0d0 - added mrw 10/May/06
        if(xx .eq. 1.0d0) then
	  ctlhpardis = 0.0d0
	  return
	endif
	
c skip the initialization in x if same as in the previous call.
	if(x .eq. xlast) goto 100
	xlast = x

      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)','Severe error: x <= 0 in ParDis x=', x 
        Stop 
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
         Msg = '0 < X < Xmin in ParDis; extrapolation used!'
         CALL CtLhWARNR (IWRN1, Msg, 'X', X, Xmin, 1D0, 1)
      Elseif (JLx .LE. Nx-2) Then
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)','Severe error: x > 1 in ParDis x=', x 
        Stop 
      Endif
      ss = x**xpow
      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)
      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4
      sy2 = ss - svec2 
      sy3 = ss - svec3 
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12 
      const6 = (s1213*sy2-s12*sy3)*tmp/s34
      EndIf

100	continue

c skip the initialization in q if same as in the previous call.
        if(q .eq. qlast) goto 110
	qlast = q

      tt = log(log(Q/Al))

      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (Q .GE. QV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif
      If     (JLq .LE. 0) Then
         Jq = 0
         If (JLq .LT. 0) Then
          Msg = 'Q < Q0 in ParDis; extrapolation used!'
          CALL CtLhWARNR (IWRN2, Msg, 'Q', Q, Qini, 1D0, 1)
         EndIf
      Elseif (JLq .LE. Nt-2) Then
         Jq = JLq - 1
      Else
        Jq = Nt - 3
        If (JLq .GE. Nt) Then
         Msg = 'Q > Qmax in ParDis; extrapolation used!'
         CALL CtLhWARNR (IWRN3, Msg, 'Q', Q, Qmax, 1D0, 1)
        Endif
      Endif

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)
      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4
      ty2 = tt - tvec2
      ty3 = tt - tvec3
      tmp1 = t12 + t13
      tmp2 = t24 + t34
      tdet = t12*t34 - tmp1*tmp2
      EndIf

110	continue

      jtmp = ((IPRTN + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1
      Do it = 1, nqvec
         J1  = jtmp + it*(NX+1) 
       If (Jx .Eq. 0) Then
         fij(1) = 0
         fij(2) = Upd(J1+1,iset) * Xa(1,2)
         fij(3) = Upd(J1+2,iset) * Xa(2,2)
         fij(4) = Upd(J1+3,iset) * Xa(3,2)
         Call CtLhPolint4 (XVpow(0), Fij(1), 4, ss, Fx, Dfx) 
         
         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 
       ElseIf  (JLx .Eq. Nx-1) Then
        Call CtLhPolint4 (XVpow(Nx-3), Upd(J1,iset), 4, ss, Fx, Dfx)
        Fvec(it) = Fx
       Else 
         sf2 = Upd(J1+1,iset)
         sf3 = Upd(J1+2,iset)
         Fvec(it) = (const5*(Upd(J1,iset) 
     &	                     - sf2*const1 + sf3*const2) 
     &               + const6*(Upd(J1+3,iset) 
     &	                     + sf2*const3 - sf3*const4) 
     &               + sf2*sy3 - sf3*sy2) / s23
       Endif
      enddo
      If (JLq .LE. 0) Then
        Call CtLhPolint4 (TV(0), Fvec(1), 4, tt, ff, Dfq)
      ElseIf (JLq .GE. Nt-1) Then
        Call CtLhPolint4 (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
        tf2 = fvec(2)
        tf3 = fvec(3)
        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23
        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 
     &	  +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)
        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf
      CtLhPARDIS = ff
      Return
      End

      SUBROUTINE CtLhPARPDF (IACT, NAME, VALUE, IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER NAME*(*), Uname*10
      LOGICAL START1
      DATA ILEVEL, LRET / 1, 1 /
      JRET = IRET
      CALL CtLhUPC (NAME, Ln, Uname)
      IF (IACT .EQ. 0 .OR. IACT .EQ. 4) then
c     >    IVALUE = NINT (VALUE)   !tentatively remove this since it seems not to be used
       print *,'Fatal error: iact=',iact
       stop
      ENDIF
      START1 = (IACT .NE. 1) .AND. (IACT .NE. 2)
c prepare to remove this stuff, since I think IACT=1 or 2 always
      if(start1) then
        print *,'Fatal error: start1=',start1
        stop
      endif
      IF (START1)  ILEVEL = 1
      GOTO (1, 2), ILEVEL
    1 START1 = .TRUE.
      ILEVEL = 0
      CALL CtLhParQcd (IACT, Uname(1:Ln), VALUE, JRET)
              IF (JRET .EQ. 1)  GOTO 11
              IF (JRET .EQ. 2)  GOTO 12
              IF (JRET .EQ. 3)  GOTO 13
              IF (JRET .GT. 4)  GOTO 15
              ILEVEL =  ILEVEL + 1
    2 CALL CtLhEVLPAR (IACT, Uname(1:Ln), VALUE, JRET)
              IF (JRET .EQ. 1)  GOTO 11
              IF (JRET .EQ. 2)  GOTO 12
              IF (JRET .EQ. 3)  GOTO 13
              IF (JRET .GT. 4)  GOTO 15
              ILEVEL =  ILEVEL + 1
      IF (.NOT. START1) GOTO 1
      IF (JRET .EQ. 0)  GOTO 10
    9 CONTINUE
      GOTO 14
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
   13 CONTINUE
   14 CONTINUE
   15 CONTINUE
      IF (JRET .NE. 4) LRET = JRET
      IF (LRET.EQ.0 .OR. LRET.EQ.2 .OR. LRET.EQ.3) THEN
        PRINT *, 'Error in CtLhPARPDF: IRET, IACT, NAME, VALUE =',
     >  LRET, IACT, NAME, VALUE
	PRINT *, 'fatal error in CtLhparpdf'
	stop
      EndIf
      IRET= JRET
      RETURN
  100 FORMAT (/)
      END
      SUBROUTINE CtLhParQcd(IACT,NAME,VALUE,IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER IACT,IRET
      CHARACTER*(*) NAME
      IRET=1
      IF (IACT.EQ.0) THEN
         WRITE (NINT(VALUE), *)  'LAM(BDA), NFL, ORD(ER), Mi, ',
     >               '(i in 1 to 9), LAMi (i in 1 to NFL)'
         IRET=4
      ELSEIF (IACT.EQ.1) THEN
         CALL CtLhQCDSET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.2) THEN
         CALL CtLhQCDGET (NAME,VALUE,IRET)
      ELSE
         IRET=3
      ENDIF
      RETURN
      END
      FUNCTION CtLhPFF1 (XX)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LA, LB, LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      PARAMETER (MX = 3)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtKRNL00 / DZ, XL(MX), NNX
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtKRN1ST / FF1(0:MXX),FG1(0:MXX),GF1(0:MXX),GG1(0:MXX),
     >                  FF2(0:MXX), FG2(0:MXX), GF2(0:MXX), GG2(0:MXX),
     >                  PNSP(0:MXX), PNSM(0:MXX)
      SAVE
      DATA LA, LB / 2 * .FALSE. /
      LB = .TRUE.
      ENTRY CtLhTFF1(ZZ)
      LA = .TRUE.
    2 IF (LA .AND. .NOT.LB) THEN
        Z = ZZ
        X = CtLhXFRMZ (Z)
      Else
        X = XX
      EndIf
      IF (X .GE. D1) THEN
        CtLhPFF1 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = CtLhZFRMX (X)
        TEM = CtLhFINTRP (FF1,  -DZ, DZ, NX,  Z,  ERR, IRT)
      Else
        CALL CtLhPOLIN1 (XL, FF1(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         IF (LB) THEN
            CtLhPFF1 = TEM / (1.-X)
            LB   =.FALSE.
         Else
            CtLhTFF1 = TEM / (1.-X) * CtLhDXDZ(Z)
         EndIf
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QFF1 = TEM
            LB   =.FALSE.
         Else
            RFF1 = TEM * X / (1.-X)
         EndIf
      EndIf
      RETURN
      ENTRY CtLhFNSP (XX)
      X = XX
      IF (X .GE. D1) THEN
        CtLhFNSP = 0.
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = CtLhZFRMX (X)
        TEM = CtLhFINTRP (PNSP,  -DZ, DZ, NX,  Z,  ERR, IRT)
      Else
        CALL CtLhPOLIN1 (XL, PNSP(1), MX, X, TEM, ERR)
      EndIf
      CtLhFNSP = TEM / (1.- X)
      RETURN
      ENTRY CtLhFNSM (XX)
      X = XX
      IF (X .GE. D1) THEN
        CtLhFNSM = 0.
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = CtLhZFRMX (X)
        TEM = CtLhFINTRP (PNSM,  -DZ, DZ, NX,  Z,  ERR, IRT)
      Else
        CALL CtLhPOLIN1 (XL, PNSM(1), MX, X, TEM, ERR)
      EndIf
      CtLhFNSM = TEM / (1.- X)
      RETURN
      ENTRY CtLhRGG1 (XX)
      X = XX
      IF (X .GE. D1) THEN
        PGG1= 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = CtLhZFRMX (X)
        TEM = CtLhFINTRP (GG1,  -DZ, DZ, NX,  Z,  ERR, IRT)
      Else
        CALL CtLhPOLIN1 (XL, GG1(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PGG1 = TEM / X / (1.-X)
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QGG1 = TEM / X
            LB   =.FALSE.
         Else
            CtLhRGG1 = TEM / (1.-X)
         EndIf
      EndIf
      RETURN
      ENTRY CtLhRFF2 (XX)
      X = XX
      IF (X .GE. D1) THEN
        PFF2 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = CtLhZFRMX (X)
        TEM = CtLhFINTRP (FF2,  -DZ, DZ, NX,  Z,  ERR, IRT)
      Else
        CALL CtLhPOLIN1 (XL, FF2(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PFF2 = TEM / X / (1.-X)
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QFF2 = TEM / X
            LB   =.FALSE.
         Else
            CtLhRFF2 = TEM / (1.-X)
         EndIf
      EndIf
      RETURN
      ENTRY CtLhRGG2 (XX)
      X = XX
      IF (X .GE. D1) THEN
        PGG2 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = CtLhZFRMX (X)
        TEM = CtLhFINTRP (GG2,  -DZ, DZ, NX,  Z,  ERR, IRT)
      Else
        CALL CtLhPOLIN1 (XL, GG2(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PGG2 = TEM / X / (1.-X)
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QGG2 = TEM / X
            LB   =.FALSE.
         Else
            CtLhRGG2 = TEM / (1.-X)
         EndIf
      EndIf
      RETURN
      END
      SUBROUTINE CtLhPOLIN1 (XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
      SUBROUTINE CtLhQARRAY (NINI)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      COMMON / LhCtQARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / LhCtQARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      NCNT = 0
      IF (NT .GE. mxq) NT = mxq - 1
      S = LOG(QINI/AL)
      TINI = LOG(S)
      S = LOG(QMAX/AL)
      TMAX = LOG(S)
    1 DT0 = (TMAX - TINI) / float(NT)
      NINI = LhCtNFL(QINI)
      NFMX = LhCtNFL(QMAX)
      Call CtLhParQcd (2, 'ORDER', Ord, Ir)
      Call CtLhParQcd (2, 'ALAM', Al0, Ir)
      Call CtLhParQcd (2, 'NFL', Afl0, Ir)
      AFL = NfMx
      Call CtLhParQcd (1, 'NFL', AFL, Ir)
      Iordr = Nint (Ord)
      Ifl0  = Nint (Afl0)
      Call CtLhSetLam (Ifl0, Al0, Iordr)
      NG = NFMX - NINI + 1
      QIN  = QINI
      QOUT = QINI
      S = LOG(QIN/AL)
      TIN  = LOG(S)
      TLN(1) = TIN
      NTL(1)  = 0
      QV(0) = QINI
      TV(0) = Tin
      DO 20 NEFF = NINI, NFMX
        ICNT = NEFF - NINI + 1
        IF (NEFF .LT. NFMX) THEN
          THRN = CtLhAMHATF (NEFF + 1)
          QOUN = MIN (QMAX, THRN)
        Else
          QOUN = QMAX
        EndIf
        IF (QOUN-QOUT .LE. 0.0001) THEN
          DT   = 0
          NITR = 0
        Else
          QOUT = QOUN
          S = LOG(QOUT/AL)
          TOUT = LOG(S)
          TEM = TOUT - TIN
          NITR = INT (TEM / DT0) + 1
          DT  = TEM / NITR
        EndIf
        DTN (ICNT) = DT
        NTN (ICNT) = NITR
        TLN (ICNT) = TIN
        NTL (ICNT+1) = NTL(ICNT) + NITR
        IF (NITR .NE. 0) THEN
        DO 205 I = 1, NITR
           TV (NTL(ICNT)+I) = TIN + DT * I
           S = EXP (TV(NTL(ICNT)+I))
           QV (NTL(ICNT)+I) = AL * EXP (S)
  205   CONTINUE
        EndIf
        QIN = QOUT
        TIN = TOUT
   20 CONTINUE
      NCNT = NCNT + 1
      NTP = NTL (NG + 1)
      ND  = NTP - NT
      IF (NTP .GE. MXQ) THEN
         NT = MXQ - ND - NCNT
         GOTO 1
      EndIf
      NT = NTP
      RETURN
      END
      SUBROUTINE CtLhQCDGET(NAME,VALUE,IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) NAME
      COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      COMMON / LhCtCOMQMS / VALQMS(9)
      LOGICAL SET
      PARAMETER (PI=3.141592653589793d0, EULER=0.57721566)
      ICODE = LhCtNAMQCD(NAME)
      IRET = 1
      IF (ICODE .EQ. 1) THEN
         VALUE = AL
      ELSEIF (ICODE .EQ. 2) THEN
         VALUE = NF
      ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 12))  THEN
         VALUE = VALQMS(ICODE - 2)
      ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
         VALUE = ALAM(ICODE - 13)
      ELSEIF (ICODE .EQ. 24) THEN
         VALUE = NORDER
      ELSE
         IRET=0
      ENDIF
      END
      SUBROUTINE CtLhQCDSET (NAME,VALUE,IRET)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) NAME
      COMMON / LhCtCOMQMS / VALQMS(9)
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      PARAMETER (PI=3.141592653589793d0, EULER=0.57721566)
      IVALUE = NINT(VALUE)
      ICODE  = LhCtNAMQCD(NAME)
      IF (ICODE .EQ. 0) THEN
         IRET=0
c     print *,'warning empty CtLhQCDSET call: NAME=',
c     &	         NAME,' VALUE=',VALUE
      ELSE
         IRET = 1
         SET = .FALSE.
         IF (ICODE .EQ. 1) THEN
            IF (VALUE.LE.0) GOTO 12
            AL=VALUE
         ELSEIF (ICODE .EQ. 2) THEN
            IF ( (IVALUE .LT. 0) .OR. (IVALUE .GT. 9)) GOTO 12
            NF = IVALUE
         ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 11))  THEN
            IF (VALUE .LT. 0) GOTO 12
            Scle = Min (Value , VALQMS(ICODE - 2))
            AlfScle = CtLhALPI(Scle) * Pi
            VALQMS(ICODE - 2) = VALUE
            Call CtLhAlfSet (Scle, AlfScle)
         ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
            IF (VALUE .LE. 0) GOTO 12
            CALL CtLhSETL1 (ICODE-13, VALUE)
         ELSEIF (ICODE .EQ. 24)  THEN
            IF ((IVALUE .LT. 1) .OR. (IVALUE .GT. 2)) GOTO 12
            NORDER = IVALUE
         ENDIF
         IF (.NOT. SET) CALL CtLhLAMCWZ
      ENDIF
      RETURN
 12   IRET=2
      RETURN
      END
      FUNCTION CtLhQZBRNT(FUNC, X1, X2, TOLIN, IRT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)
      external func
      TOL = ABS(TOLIN)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        WRITE (*, *) 'Root must be bracketed for CtLhQZBRNT.'
        IRT = 1
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          CtLhQZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE (*, *) 'CtLhQZBRNT exceeding maximum iterations.'
      IRT = 2
      CtLhQZBRNT=B
      RETURN
      END
      SUBROUTINE CtLhRATINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
      FUNCTION CtLhRTALF (EFLLN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      include 'parmsetup.inc'
      PARAMETER (PI = 3.141592653589793d0)
      COMMON / CtLhRTALFC / ALFST, JORD, NEFF
      EFMULM = EXP (EFLLN)
      TEM1 = PI / ALFST
      TEM2 = 1. / CtLhALPQCD (JORD, NEFF, EFMULM, I)
      CtLhRTALF = TEM1 - TEM2
      END
      Subroutine CtLhbldat1
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      include 'parmsetup.inc'
      LOGICAL LSTX
      PARAMETER (MXX = 105, MXQ = 25, MxF = 6)
      PARAMETER (MxPN = MxF * 2 + 2)
      PARAMETER (MxQX= MXQ * MXX,   MxPQX = MxQX * MxPN)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtQARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / LhCtPEVLDT / UPD(MXPQX,nmxset), KF, Nelmt
	PARAMETER (NF0 = 4, Nshp = 8,NEX = Nshp+2)
      XMIN =  .999999D-4
      XCR = 1.5 
      JT = 1
      Return
      END
      SUBROUTINE CtLhSETL1  (NEF, VLAM)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL SET
      COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      COMMON / LhCtCOMQMS / VALQMS(9)
      IF (NEF .LT. 0 .OR. NEF .GT. NF) THEN
        WRITE(*,*)'NEF out of range in CtLhSETL1: NEF NF =',NEF,NF
        STOP
      ENDIF
      AMHAT(0) = 0.
      DO 5 N = 1, NF
         AMHAT(N) = VALQMS(N)
    5    CONTINUE
      ALAM(NEF) = VLAM
      DO 10 N = NEF, 1, -1
         CALL CtLhTRNLAM(NORDER, N, -1, IR1)
   10    CONTINUE
      DO 20 N = NEF, NF-1
         CALL CtLhTRNLAM(NORDER, N, 1, IR1)
   20    CONTINUE
      DO 30, N = NF, 1, -1
         IF ((ALAM(N) .GE. 0.7 * AMHAT(N))
     >       .OR. (ALAM(N-1) .GE. 0.7 * AMHAT(N)))THEN
            NHQ = NF - N
            GOTO 40
            ENDIF
   30    CONTINUE
      NHQ = NF
   40 CONTINUE
      DO 50, N = NF-NHQ, 1, -1
         AMHAT(N) = 0
         ALAM(N-1) = ALAM(N)
   50    CONTINUE
      AMN = ALAM(NF)
      DO 60, N = 0, NF-1
         IF (ALAM(N) .GT. AMN)  AMN = ALAM(N)
   60    CONTINUE
      AMN = AMN * 1.0001
      AL = ALAM(NF)
      SET = .TRUE.
      RETURN
      END
      SUBROUTINE CtLhSETLAM (NEF, WLAM, IRDR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      IF ((NEF .LT. 0) .OR. (NEF .GT. NF)) THEN
         WRITE(*,*)'NEF out of range in CtLhSETLAM: NEF NF=',NEF,NF
         STOP
      ENDIF
      VLAM = WLAM
      IF (IRDR .NE. NORDER) then
	PRINT *,'fatal error: wanted cnvl1'
	stop
      ENDIF
      CALL CtLhSETL1 (NEF, VLAM)
      END
      Subroutine CtLhbldat2
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtCOMQMS / VALQMS(9)
      COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
      LOGICAL SET
      AL = .226d0
      NF = 5
      NORDER = 2
      SET = .FALSE.
      VALQMS(1) =  0.
      VALQMS(2) =  0.
      VALQMS(3) =  0.2d0
      VALQMS(4) =  1.3d0
      VALQMS(5) =  4.5d0
      VALQMS(6) =  174.d0
      VALQMS(7) =  0.
      VALQMS(8) =  0.
      VALQMS(9) =  0.
      Return
      END
      FUNCTION CtLhSMPNOL (NX, DX, FN, ERR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)
      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN CtLhSMPNOL!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(2)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(2) * 2.
      ELSE
         IF (MS .EQ. 0) THEN
            TEM = DX * (23.* FN(2) - 16.* FN(3) + 5.* FN(4)) / 12.
            TMP = DX * (3.* FN(2) - FN(3)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + CtLhSMPSNA (NX-1, DX, FN(2), ER1)
            ERR = ABS(ER1) + ERR
         ELSE
            TEM = DX * (8.* FN(2) - 4.* FN(3) + 8.* FN(4)) / 3.
            TMP = DX * (3.* FN(2) + 2.* FN(3) + 3.* FN(4)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + CtLhSMPSNA (NX-4, DX, FN(5), ER1)
            ERR = ABS(ER1) + ERR
         ENDIF
      ENDIF
      CtLhSMPNOL = TEM
      RETURN
      END
      FUNCTION CtLhSMPSNA (NX, DX, F, ERR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MAXX = 1000)
      DIMENSION F(NX)
      DATA IW1, IW2, TINY / 2*0, 1.E-35 /
      IF (DX .LE. 0.) THEN
        CALL CtLhWARNR(IW2,'DX cannot be < 0. in CtLhSMPSNA', 'DX', 
     >         DX, D0, D1, 0)
        CtLhSMPSNA = 0.
        RETURN
      ENDIF
      IF (NX .LE. 0 .OR. NX .GT. MAXX) THEN
        CALL CtLhWARNI(IW1, 'NX out of range in CtLhSMPSNA', 'NX', NX,
     >               1, MAXX, 1)
        SIMP = 0.
      ELSEIF (NX .EQ. 1) THEN
        SIMP = 0.
      ELSEIF (NX .EQ. 2) THEN
        SIMP = (F(1) + F(2)) / 2.
        ERRD = (F(1) - F(2)) / 2.
      ELSE
        MS = MOD(NX, 2)
        IF (MS .EQ. 0) THEN
          ADD = (9.*F(NX) + 19.*F(NX-1) - 5.*F(NX-2) + F(NX-3)) / 24.
          NZ = NX - 1
        ELSE
          ADD = 0.
          NZ = NX
        ENDIF
        IF (NZ .EQ. 3) THEN
          SIMP = (F(1) + 4.* F(2) + F(3)) / 3.
          TRPZ = (F(1) + 2.* F(2) + F(3)) / 2.
        ELSE
          SE = F(2)
          SO = 0
          NM1 = NZ - 1
          DO 60 I = 4, NM1, 2
            IM1 = I - 1
            SE = SE + F(I)
            SO = SO + F(IM1)
   60     CONTINUE
          SIMP = (F(1) + 4.* SE + 2.* SO + F(NZ)) / 3.
          TRPZ = (F(1) + 2.* (SE + SO) + F(NZ)) / 2.
        ENDIF
        ERRD = TRPZ - SIMP 
        SIMP = SIMP + ADD
      ENDIF
      CtLhSMPSNA = SIMP * DX
      IF (ABS(SIMP) .GT. TINY) THEN
        ERR = ERRD / SIMP
      ELSE
        ERR = 0.
      ENDIF
      RETURN
      END
      SUBROUTINE CtLhSNEVL(IKNL,NX,NT,JT,DT,TIN,NEFF,UI,GI,US,GS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXQX= MXQ * MXX)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      DIMENSION UI(NX), US(0:NX, 0:NT)
      DIMENSION GI(NX), GS(0:NX, 0:NT)
      DIMENSION Y0(MXX), Y1(MXX), YP(MXX), F0(MXX), F1(MXX), FP(MXX)
      DIMENSION Z0(MXX), Z1(MXX), ZP(MXX), G0(MXX), G1(MXX), GP(MXX)
      DATA D0 / 0.0 /
      JTT = 2 * JT
      DDT = DT / JTT
      IF (NX .GT. MXX) THEN
      WRITE (*,*) 'Nx =', NX, ' too many pts in CtLhSNEVL'
      STOP 'Program stopped in CtLhSNEVL'
      EndIf
      TMD = TIN + DT * NT / 2.
      AMU = EXP(TMD)
      TEM = 6./ (33.- 2.* NEFF) / CtLhALPI(AMU)
      TLAM = TMD - TEM
      DO 9 IX = 1, NX
      US (IX, 0) = UI(IX)
      GS (IX, 0) = GI(IX)
    9 CONTINUE
      US ( 0, 0) = (UI(1) - UI(2))* 3D0 + UI(3)
      GS ( 0, 0) = (GI(1) - GI(2))* 3D0 + GI(3)
      TT = TIN
      DO 10 IZ = 1, NX
      Y0(IZ) = UI(IZ)
      Z0(IZ) = GI(IZ)
   10 CONTINUE
      DO 20 IS = 1, NT
         DO 202 JS = 1, JTT
            IRND = (IS-1) * JTT + JS
            IF (IRND .EQ. 1) THEN
                CALL CtLhSNRHS (TT, NEFF, Y0,Z0,  F0,G0)
                DO 250 IZ = 1, NX
                   Y0(IZ) = Y0(IZ) + DDT * F0(IZ)
                   Z0(IZ) = Z0(IZ) + DDT * G0(IZ)
  250           CONTINUE
                TT = TT + DDT
                CALL CtLhSNRHS (TT, NEFF, Y0, Z0,  F1, G1)
                DO 251 IZ = 1, NX
                   Y1(IZ) = UI(IZ) + DDT * (F0(IZ) + F1(IZ)) / 2D0
                   Z1(IZ) = GI(IZ) + DDT * (G0(IZ) + G1(IZ)) / 2D0
  251           CONTINUE
            Else
                CALL CtLhSNRHS (TT, NEFF, Y1, Z1,  F1, G1)
                DO 252 IZ = 1, NX
                   YP(IZ) = Y1(IZ) + DDT * (3D0 * F1(IZ) - F0(IZ)) / 2D0
                   ZP(IZ) = Z1(IZ) + DDT * (3D0 * G1(IZ) - G0(IZ)) / 2D0
  252           CONTINUE
                TT = TT + DDT
                CALL CtLhSNRHS (TT, NEFF, YP, ZP,  FP, GP)
                DO 253 IZ = 1, NX
                   Y1(IZ) = Y1(IZ) + DDT * (FP(IZ) + F1(IZ)) / 2D0
                   Z1(IZ) = Z1(IZ) + DDT * (GP(IZ) + G1(IZ)) / 2D0
                   F0(IZ) = F1(IZ)
                   G0(IZ) = G1(IZ)
  253           CONTINUE
            EndIf
  202    CONTINUE
         DO 260 IX = 1, NX
           IF (IKNL .GT. 0) THEN
            US (IX, IS) = MAX(Y1(IX), D0)
            GS (IX, IS) = MAX(Z1(IX), D0)
           Else
            US (IX, IS) = Y1(IX)
            GS (IX, IS) = Z1(IX)
           EndIf
  260    CONTINUE
         US(0, IS) = 3D0*Y1(1) - 3D0*Y1(2) + Y1(3)
         GS(0, IS) = 3D0*Z1(1) - 3D0*Z1(2) + Z1(3)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CtLhSNRHS (TT, NEFF, FI, GI,  FO, GO)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtXYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / LhCtKRNL01 / AFF2(MXX),AFG2(MXX),AGF2(MXX),AGG2(MXX),
     >                  ANSP (MXX), ANSM (MXX), ZFG2, ZGF2, ZQQB
      COMMON / LhCtKRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / LhCtEVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      DIMENSION GI(NX), GO(NX), G1(MXX), G2(MXX), G3(MXX), G4(MXX)
      DIMENSION FI(NX), FO(NX), W0(MXX), W1(MXX), WH(MXX), WM(MXX)
      DIMENSION R0(MXX), R1(MXX), R2(MXX), RH(MXX), RM(MXX)
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = CtLhALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
      CALL CtLhINTEGR (NX,-1, FI, WM, IR1)
      CALL CtLhINTEGR (NX, 0, FI, W0, IR2)
      CALL CtLhINTEGR (NX, 1, FI, W1, IR3)
      CALL CtLhINTEGR (NX,-1, GI, RM, IR4)
      CALL CtLhINTEGR (NX, 0, GI, R0, IR5)
      CALL CtLhINTEGR (NX, 1, GI, R1, IR6)
      CALL CtLhINTEGR (NX, 2, GI, R2, IR7)
      CALL CtLhHINTEG (NX,    FI, WH)
      CALL CtLhHINTEG (NX,    GI, RH)
      IF (IKNL .GT. 0) THEN
      DO 230 IZ = 1, NX
      FO(IZ) = ( 2D0 * FI(IZ)
     >      + 4D0 / 3D0 * ( 2D0 * WH(IZ) - W0(IZ) - W1(IZ) ))
     >      + NEFF * ( R0(IZ) - 2D0 * R1(IZ) + 2D0 * R2(IZ) )
      FO(IZ) = FO(IZ) * CPL
      GO(IZ) = 4D0 / 3D0 * ( 2D0 * WM(IZ) - 2D0 * W0(IZ)  + W1(IZ) )
     >      + (33D0 - 2D0 * NEFF) / 6D0 * GI(IZ)
     >      + 6D0 * (RH(IZ) + RM(IZ) - 2D0 * R0(IZ) + R1(IZ) - R2(IZ))
      GO(IZ) = GO(IZ) * CPL
  230 CONTINUE
      Else
      DO 240 IZ = 1, NX
      FO(IZ) = NEFF * (-R0(IZ) + 2.* R1(IZ) )
     > + 2.* FI(IZ) + 4./ 3.* ( 2.* WH(IZ) - W0(IZ) - W1(IZ) )
      FO(IZ) = FO(IZ) * CPL
      GO(IZ) = 4./ 3.* ( 2.* W0(IZ) - W1(IZ) )
     >+ (33.- 2.* NEFF) / 6.* GI(IZ) + 6.*(RH(IZ) + R0(IZ) - 2.* R1(IZ))
      GO(IZ) = GO(IZ) * CPL
  240 CONTINUE
      EndIf
      IF (IKNL .EQ. 2) THEN
      DZ = 1./(NX - 1)
      DO 21 I = 1, NX-1
        X = XV(I)
        NP = NX - I + 1
        IS = NP
           g2(1)=0d0
           g3(1)=0d0
        DO 31 KZ = 2, NP
          IY = I + KZ - 1
          IT = NX - IY + 1
          XY = ZZ (IS, IT)
          G1(KZ) = FFG(I, IY) * (FI(IY) - XY**2 *FI(I))
          G4(KZ) = GGF(I, IY) * (GI(IY) - XY**2 *GI(I))
           G2(KZ) = FFG(IS,IT) * (GI(IY) - xy*GI(I))    !FG
           G3(KZ) = GGF(IS,IT) * (FI(IY) - XY*FI(I))    !GF (usual notations)
   31   CONTINUE
        TEM1 = CtLhSMPNOL (NP, DZ, G1, ERR)
        TEM2 = CtLhSMPSNA (NP, DZ, G2, ERR)
        TEM3 = CtLhSMPSNA (NP, DZ, G3, ERR)
        TEM4 = CtLhSMPNOL (NP, DZ, G4, ERR)
        TEM1 = TEM1 - FI(I) * (AFF2(I) + ZGF2)
        TEM4 = TEM4 - GI(I) * (AGG2(I) + ZFG2)
         tem2 = tem2 + GI(I)*AFG2(I)
         tem3=  tem3 + FI(I)*AGF2(I)
        TMF = TEM1 + TEM2
        TMG = TEM3 + TEM4
        FO(I) = FO(I) + TMF * CPL2
        GO(I) = GO(I) + TMG * CPL2
   21 CONTINUE
      EndIf
      RETURN
      END
      FUNCTION CtLhSPENC2 (X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL CtLhSPN2IN
      COMMON / LhCtSPENCC / XX
      DATA U1, AERR, RERR / 1.D0, 1.E-7, 5.E-3 /
      XX = X
      TEM = CtLhGausInt(CtLhSPN2IN, XX, U1, AERR, RERR, ERR, IRT)
      CtLhSPENC2 = TEM + LOG (XX) ** 2 / 2.
      RETURN
      END
      FUNCTION CtLhSPN2IN (ZZ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtSPENCC / X
      Z = ZZ
      TEM = LOG (1.+ X - Z) / Z
      CtLhSPN2IN = TEM
      RETURN
      END
      SUBROUTINE CtLhSTUPKL (NFL)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MX = 3)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtXYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtKRN1ST / FF1(0:MXX),FG1(0:MXX),GF1(0:MXX),GG1(0:MXX),
     >                  FF2(0:MXX), FG2(0:MXX), GF2(0:MXX), GG2(0:MXX),
     >                  PNSP(0:MXX), PNSM(0:MXX)
      COMMON / LhCtKRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / LhCtKRNL00 / DZ, XL(MX), NNX
      COMMON / LhCtKRNL01 / AFF2(MXX),AFG2(MXX),AGF2(MXX),AGG2(MXX),
     >                  ANSP (MXX), ANSM (MXX), ZFG2, ZGF2, ZQQB
      EXTERNAL CtLhPFF1, CtLhRGG1, CtLhRFF2, CtLhRGG2
      EXTERNAL CtLhFNSP, CtLhFNSM
      dimension aff1(mxx),agg1(mxx)
      PARAMETER (PI = 3.141592653589793d0, PI2 = PI**2)
      DATA CF, CG, TR / 1.33333333333333d0, 3.0, 0.5 /
      data zeta3/1.20205690315959d0/          ! zeta(3.0)
      SAVE
      DATA AERR, RERR / 0.0, 0.02 /
      NNX = NX
      DZ = 1./ (NX - 1)
      DO 5 I0 = 1, MX
        XL(I0) = XV(I0)
    5 CONTINUE
      DO 10 I = 1, NX-1
        XZ = XV(I)
      CALL CtLhKERNEL (XZ, FF1(I), GF1(I), FG1(I), GG1(I), PNSP(I),
     >          PNSM(I), FF2(I), GF2(I), FG2(I), GG2(I), NFL, IRT)
   10 CONTINUE
      FF1(0) = FF1(1) * 3. - FF1(2) * 3. + FF1(3)
      FG1(0) = FG1(1) * 3. - FG1(2) * 3. + FG1(3)
      GF1(0) = GF1(1) * 3. - GF1(2) * 3. + GF1(3)
      GG1(0) = GG1(1) * 3. - GG1(2) * 3. + GG1(3)
      PNSP(0) = PNSP(1) * 3. - PNSP(2) * 3. + PNSP(3)
      PNSM(0) = PNSM(1) * 3. - PNSM(2) * 3. + PNSM(3)
      FF2(0) = FF2(1) * 3. - FF2(2) * 3. + FF2(3)
      FG2(0) = FG2(1) * 3. - FG2(2) * 3. + FG2(3)
      GF2(0) = GF2(1) * 3. - GF2(2) * 3. + GF2(3)
      GG2(0) = GG2(1) * 3. - GG2(2) * 3. + GG2(3)
      FF1(NX) = FF1(NX-1) * 3. - FF1(NX-2) * 3. + FF1(NX-3)
      FG1(NX) = FG1(NX-1) * 3. - FG1(NX-2) * 3. + FG1(NX-3)
      GF1(NX) = GF1(NX-1) * 3. - GF1(NX-2) * 3. + GF1(NX-3)
      GG1(NX) = GG1(NX-1) * 3. - GG1(NX-2) * 3. + GG1(NX-3)
      PNSM(NX) = PNSM(NX-1) * 3. - PNSM(NX-2) * 3. + PNSM(NX-3)
      PNSP(NX) = PNSP(NX-1) * 3. - PNSP(NX-2) * 3. + PNSP(NX-3)
      FF2(NX) = FF2(NX-1) * 3. - FF2(NX-2) * 3. + FF2(NX-3)
      FG2(NX) = FG2(NX-1) * 3. - FG2(NX-2) * 3. + FG2(NX-3)
      GF2(NX) = GF2(NX-1) * 3. - GF2(NX-2) * 3. + GF2(NX-3)
      GG2(NX) = GG2(NX-1) * 3. - GG2(NX-2) * 3. + GG2(NX-3)
         RER = RERR * 4.
         AFF1(1) = CtLhGausInt(CtLhPFF1,D0,XV(1),AERR,RERR,ER1,IRT)
         DGG1     = NFL / 3.
         TMPG     = CtLhGausInt(CtLhRGG1,D0,XV(1),AERR,RERR,ER3,IRT)
         AGG1(1) = TMPG + DGG1
       ANSM(1) = CtLhGausInt(CtLhFNSM,D0,XV(1),AERR,RER,ER2,IRT)
       ANSP(1) = CtLhGausInt(CtLhFNSP,D0,XV(1),AERR,RER,ER2,IRT)
         AER = AFF1(1) * RER
         AFF2(1) = CtLhGausInt(CtLhRFF2, D0, XV(1),  AER, RER, ER2, IRT)
         AER = AGG1(1) * RER
         AGG2(1) = CtLhGausInt(CtLhRGG2, D0, XV(1),  AER, RER, ER4, IRT)
      DO 20 I2 = 2, NX-1
      TEM =CtLhGausInt(CtLhPFF1,XV(I2-1),XV(I2),AERR,RERR,ER1,IRT)
      AFF1(I2) = TEM + AFF1(I2-1)
      AER = ABS(TEM * RER)
      AFF2(I2)=CtLhGausInt(CtLhRFF2,XV(I2-1),XV(I2),AER,RER,ER2,IRT)
     >        +AFF2(I2-1)
      TEM      = CtLhGausInt(CtLhRGG1,XV(I2-1),XV(I2),AERR,RERR,ER3,IRT)
      TMPG     = TMPG + TEM
      AGG1(I2) = TMPG + DGG1
      AER = ABS(TEM * RER)
      AGG2(I2)=CtLhGausInt(CtLhRGG2,XV(I2-1),XV(I2),AER,RER,ER4,IRT)
     >        +AGG2(I2-1)
      ANSP(I2)=CtLhGausInt(CtLhFNSP,XV(I2-1),XV(I2),AERR,RER,ER4,IRT)
     >        +ANSP(I2-1)
      ANSM(I2)=CtLhGausInt(CtLhFNSM,XV(I2-1),XV(I2),AERR,RER,ER4,IRT)
     >        +ANSM(I2-1)
   20 CONTINUE
      ANSP(NX)=CtLhGausInt(CtLhFNSP,XV(NX-1),D1,AERR,RER,ERR,
     > IRT) + ANSP(NX-1)
      ANSM(NX)=CtLhGausInt(CtLhFNSM,XV(NX-1),D1,AERR,RER,ERR,
     > IRT) + ANSM(NX-1)
           TRNF = TR * NFL
      do i2=1,nx-1                                        !loop over x
         x=xv(i2)
         XI = 1./ X                                    !auxiliary definitions
         X2 = X ** 2
         X3=  x**3
         XLN = DLOG (X)
         XLN2 = XLN ** 2
         XLN1M = DLOG (1.- X)
         xLi2m=CtLhxLi(2,-x)
         xLi2=CtLhxLi(2,x)
         xLi3=CtLhxLi(3,x)
         xLi31m=CtLhxLi(3,1d0-x)
         xLi32=CtLhxLi(3,x2)
         xln1m2=xln1m*xln1m
         xln1p=dlog(1d0+x)
         x1m=1d0-x
         x1p=1d0+x
         x3m=3d0-x
         x3p=3d0+x
         wgfcft=
     > (9 + 4*Pi2 - 22*x + 13*x2 + 6*(3 - 4*x + x2)*xln1m +
     > 40*xln - 24*xLi2)/9.
       wgfcf2=
     > (6*(2*(-9 + Pi2) + 3*x*(5 + x)) +4*(3 +2*Pi2+3*x*(-3 + 2*x))*
     > xln1m + 6*x3m*x1m*xln1m2 - 6*(x*(8 + 3*x) + 4*xln1m2)*
     > xln - 3*(-4 + x)*x*xln2)/12 - 2*(3 + 2*xln1m)*xLi2 - 4*xLi31m
       wgfcfg=
     > (3637-186*Pi2-x*(3198+72*Pi2+x*(231 + 208*x)))/108.- xln +
     > (3*xln1m*(-33 - 4*Pi2 + (50 - 17*x)*x - 3*x3m*x1m*xln1m) +
     > 2*(x*(198 + x*(27+8*x))+9*xln1m*(3 - 4*x + x2 + 2*xln1m))*
     > xln - 9*x*(4 + x)*xln2)/18- x1p*x3p*xln*xln1p-
     > (x1p*x3p - 4*xln)*xLi2m + (31d0/3d0 +4*xln1m- 4*xln)*xLi2 +
     > 4*xLi31m + 12*xLi3 - 2*xLi32 - 10*zeta3
       wfgcft=
     > (18 - 81*x + 6*Pi2*x + 123*x2 - 6*Pi2*x2 - 60*x3 +
     > 4*Pi2*x3 - 6*(-2 + 3*x - 3*x2 + 2*x3)*xln1m2 -33*x*xln +
     > 15*x2*xln - 24*x3*xln - 9*x*xln2 + 9*x2*xln2 -
     > 12*x3*xln2 - 12*x1m*xln1m*(-1 + 2*x2 + 2*xln - x*xln +
     > 2*x2*xln) - 24*xLi2)/9.
       wfgcgt=
     > (2*(-67 + 2*Pi2 + x*(64 + x*(-91 + 3*Pi2 + 94*x)) +
     > x1m*(7+x*(-5+16*x))*xln1m -3*x1m*(2+ x*(-1+2*x))*xln1m2 -
     > 20*xln - 3*x*xln*(13 + 16*x*x1p - 3*x1p*xln) +
     > 6*x1p*(2+x+2*x2)*xln*xln1p+6*x1p*(2+x+2*x2)*xLi2m))/9.
       AGF2(I2) = CF*TRNF*WGFCFT+CF**2* WGFCF2+CF*CG*WGFCFG
       AFG2(I2) = CF*TRNF*WFGCFT            +CG*TRNF*WFGCGT
      enddo !i2
       AGF2(nx)=0d0
       AFG2(nx)=0d0
       ZGF2=-28./27.*Cf**2+94./27.*Cf*Cg -52./27.*Cf*TrNf
       ZFG2= 37./27.*Cf*TrNf + 35./54.*Cg*TrNf
       ZQQB=1.43862321154902*(Cf**2-0.5*Cf*Cg)
      DO 21 IX = 1, NX-1
        X = XV(IX)
        NP = NX - IX + 1
        IS = NP
        XG2 = (LOG(1./(1.-X)) + 1.) ** 2
        FFG (IS, IS) = FG2(NX) * DXTZ(I) * XG2
      GGF (IS, IS) = GF2(NX) * DXTZ(I) * XG2
      PNS (IS, IS) =PNSM(NX) * DXTZ(I)
        DO 31 KZ = 2, NP
          IY = IX + KZ - 1
          IT = NX - IY + 1
          XY = X / XV(IY)
          XM1 = 1.- XY
          XG2 = (LOG(1./XM1) + 1.) ** 2
          Z  = ZZ (IX, IY)
          TZ = (Z + DZ) / DZ
          IZ = TZ
          IZ = MAX (IZ, 0)
          IZ = MIN (IZ, NX-1)
          DT = TZ - IZ
          TEM = (FF2(IZ) * (1.- DT) + FF2(IZ+1) * DT) / XM1 / XY
          FFG (IX, IY) = TEM * DXTZ(IY)
          TEM = (FG2(IZ) * (1.- DT) + FG2(IZ+1) * DT) * XG2 / XY
          FFG (IS, IT) = TEM * DXTZ(IY)
          TEM = (GF2(IZ) * (1.- DT) + GF2(IZ+1) * DT) * XG2 / XY
        GGF (IS, IT) = TEM * DXTZ(IY)
          TEM = (GG2(IZ) * (1.- DT) + GG2(IZ+1) * DT) / XM1 / XY
        GGF (IX, IY) = TEM * DXTZ(IY)
        TEM = (PNSP(IZ) * (1.- DT) + PNSP(IZ+1) * DT) / XM1
        PNS (IX, IY) = TEM * DXTZ(IY)
        TEM = (PNSM(IZ) * (1.- DT) + PNSM(IZ+1) * DT) / XM1
        PNS (IS, IT) = TEM * DXTZ(IY)
   31   CONTINUE
   21 CONTINUE
      RETURN
      END
      SUBROUTINE CtLhTRNLAM (IRDR, NF, IACT, IRT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / LhCtTRNCOM / VMULM, JRDR, N, N1
      EXTERNAL CtLhZBRLAM
      DATA ALM0, BLM0, RERR / 0.01, 10.0, 0.0001 /
      DATA IR1, SML / 0, 1.E-5 /
      IRT = 0
      N = NF
      JRDR = IRDR
      JACT = IACT
      VLAM = ALAM(N)
      IF (JACT .GT. 0) THEN
         N1 = N + 1
         THMS = AMHAT(N1)
         ALM = LOG (THMS/VLAM)
         BLM = BLM0
      ELSE
         N1 = N -1
         THMS = AMHAT(N)
         ALM = ALM0
         THMS = MAX (THMS, SML)
         BLM = LOG (THMS/VLAM)
      ENDIF
      IF (VLAM .GE. 0.7 * THMS) THEN
         IF (JACT . EQ. 1) THEN
            AMHAT(N1) = 0
         ELSE
            AMHAT(N) = 0
         ENDIF
         IRT = 4
         ALAM(N1) = VLAM
         RETURN
      ENDIF
      IF (ALM .GE. BLM) THEN
         WRITE (*, *) 'CtLhTRNLAM has ALM >= BLM: ', ALM, BLM
         WRITE (*, *) 'I do not know how to continue'
         STOP
         ENDIF
      VMULM = THMS/VLAM
      ERR = RERR * LOG (VMULM)
      WLLN = CtLhQZBRNT (CtLhZBRLAM, ALM, BLM, ERR, IR1)
      ALAM(N1) = THMS / EXP (WLLN)
      IF (IR1 .NE. 0) THEN
         WRITE (*, *) 'CtLhQZBRNT failed in CtLhTRNLAM; ',
     >        'NF, VLAM =', NF, VLAM
         WRITE (*, *) 'I do not know how to continue'
        STOP
      ENDIF
      RETURN
      END
      SUBROUTINE CtLhUPC (A, La, UpA)
      CHARACTER A*(*), UpA*(*), C*(1)
      INTEGER I, La, Ld
      La = Len(A)
      Lb = Len(UpA)
      If (Lb .Lt. La) Stop 'UpCase conversion length mismatch!'
      Ld = ICHAR('A')-ICHAR('a')
      DO 1 I = 1, Lb
        If (I .Le. La) Then
         c = A(I:I)
         IF ( LGE(C, 'a') .AND. LLE(C, 'z') ) THEN
           UpA (I:I) = CHAR(Ichar(c) + ld)
         Else
           UpA (I:I) = C
         ENDIF
        Else
         UpA (I:I) = ' '
        Endif
 1    CONTINUE
      
      RETURN
      END
      SUBROUTINE CtLhWARNI (IWRN, MSG, NMVAR, IVAB,
     >                  IMIN, IMAX, IACT)
      CHARACTER*(*) MSG, NMVAR
      Save Iw
      Data Nmax / 100 /
      IW = IWRN
      IV = IVAB
      
      IF  (IW .EQ. 0) THEN
         PRINT '(1X,A/1X, 2A,I10 /A,I4)', MSG, NMVAR, ' = ', IV
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2I10)', ' The limits are: ', IMIN, IMAX
         ENDIF
      ENDIF
      If (Iw .LT. Nmax) Then
         PRINT '(1X,A/1X,I10,A, I10)', MSG, NMVAR, ' = ', IV
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', 'CtLhWARNI Severe Warning: Too many errors'
      Endif
      IWRN = IW + 1
      RETURN
      END
      SUBROUTINE CtLhWARNR (IWRN, MSG, NMVAR, VARIAB,
     >                  VMIN, VMAX, IACT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      CHARACTER*(*) MSG, NMVAR
      Save Iw
      Data Nmax / 100 /
      IW = IWRN
      VR = VARIAB
      IF  (IW .EQ. 0) THEN
         PRINT '(1X, A/1X,2A,1PD16.7/A,I4)', MSG, NMVAR, ' = ', VR
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2(1PE15.4))', ' The limits are: ', VMIN, VMAX
         ENDIF
      ENDIF
      If (Iw .LT. Nmax) Then
         PRINT '(I5, 2A/1X,2A,I10,1PD16.7)', IW, '   ', MSG,
     >                  NMVAR, ' = ', VR
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', 'CtLhWARNR Severe Warning: Too many errors'
      Endif
      IWRN = IW + 1
      RETURN
      END
      SUBROUTINE CtLhXARRAY
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0 = 0.0, D10=10.0)
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      Character Msg*80
      COMMON / LhCtVARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / LhCtVARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
      COMMON / LhCtHINTEC / GH(NDG, MXX)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtXYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      DIMENSION G1(NDG,NDH), G2(NDG,NDH), A(NDG)
      DATA F12, F22, F32 / 1D0, 1D0, 1D0 /
      DATA (G1(I,NDH), G2(I,1), I=1,NDG) / 0.0,0.0,0.0,0.0,0.0,0.0 /
      DATA PUNY / 1D-30 /
      XV(0) = 0D0
      DZ = 1D0 / (NX-1)
      DO 10 I = 1, NX - 1
         Z = DZ * (I-1)
         ZV(I) = Z
         X = CtLhXFRMZ (Z)
         DXTZ(I) = CtLhDXDZ(Z) / X
         XV (I)  = X
         XA(I, 1) = X
         XA(I, 0) = LOG (X)
         DO 20 L = L1, L2
          IF (L .NE. 0 .AND. L .NE. 1)  XA(I, L) = X ** L
   20    CONTINUE
   10 CONTINUE
         XV(1) = Xmin
         XV(NX) = 1D0
         ZV(Nx) = 1D0
         DXTZ(NX) = CtLhDXDZ(1.D0)
         DO 21 L = L1, L2
            XA (NX, L) = 1D0
   21    CONTINUE
         XA (NX, 0) = 0D0
      DO 11 I = 1, NX-1
         ELY(I) = LOG(1D0 - XV(I))
   11 CONTINUE
       ELY(NX) = 3D0* ELY(NX-1) - 3D0* ELY(NX-2) + ELY(NX-3)
      DO 17 IX = 1, NX
      ZZ (IX, IX) = 1.
      DO 17 IY = IX+1, NX
         XY = XV(IX) / XV(IY)
         ZZ (IX, IY) = CtLhZFRMX (XY)
         ZZ (NX-IX+1, NX-IY+1) = XY
   17 CONTINUE
      DO 30 I = 1, NX-1
      IF (I .NE. NX-1) THEN
        F11 = 1D0/XV(I)
        F21 = 1D0/XV(I+1)
        F31 = 1D0/XV(I+2)
        F13 = XV(I)
        F23 = XV(I+1)
        F33 = XV(I+2)
        DET = F11*F22*F33 + F21*F32*F13 + F31*F12*F23
     >      - F31*F22*F13 - F21*F12*F33 - F11*F32*F23
        IF (ABS(DET) .LT. PUNY) THEN
           Msg='Determinant close to zero; will be arbitrarily set to:'
           CALL CtLhWARNR(IWRN, Msg, 'DET', PUNY, D0, D0, 0)
           DET = PUNY
        EndIf
        G2(1,2) = (F22*F33 - F23*F32) / DET
        G2(1,3) = (F32*F13 - F33*F12) / DET
        G2(1,4) = (F12*F23 - F13*F22) / DET
        G2(2,2) = (F23*F31 - F21*F33) / DET
        G2(2,3) = (F33*F11 - F31*F13) / DET
        G2(2,4) = (F13*F21 - F11*F23) / DET
        G2(3,2) = (F21*F32 - F22*F31) / DET
        G2(3,3) = (F31*F12 - F32*F11) / DET
        G2(3,4) = (F11*F22 - F12*F21) / DET
        B2 = LOG (XV(I+2)/XV(I))
        B3 = XV(I) * (B2 - 1.) + XV(I+2)
        GH (1,I) = B2 * G2 (2,2) + B3 * G2 (3,2)
        GH (2,I) = B2 * G2 (2,3) + B3 * G2 (3,3)
        GH (3,I) = B2 * G2 (2,4) + B3 * G2 (3,4)
      EndIf
        DO 51 J = 1, NDH
           DO 52 L = 1, NDG
              IF     (I .EQ. 1) THEN
                 GB(L,J,I) = G2(L,J)
              ElseIF (I .EQ. NX-1) THEN
                 GB(L,J,I) = G1(L,J)
              Else
                 GB(L,J,I) = (G1(L,J) + G2(L,J)) / 2D0
              EndIf
   52      CONTINUE
   51   CONTINUE
        DO 35 MM = M1, M2
           DO 40 K = 1, NDG
             KK = K + MM - 2
             IF (KK .EQ. 0) THEN
               A(K) = XA(I+1, 0) - XA(I, 0)
             Else
               A(K) = (XA(I+1, KK) - XA(I, KK)) / DBLE(KK)
             EndIf
   40      CONTINUE
           DO 41 J = 1, NDH
             TEM = 0
             DO 43 L = 1, NDG
               TEM = TEM + A(L) * GB(L,J,I)
   43        CONTINUE
             H(J,I,MM) = TEM
   41      CONTINUE
   35   CONTINUE
      DO 42 J = 1, NDG
        DO 44 L = 1, NDG
           G1(L,J) = G2(L,J+1)
   44 CONTINUE
   42 CONTINUE
   30 CONTINUE
      LSTX = .TRUE.
      RETURN
      END
      FUNCTION CtLhXFRMZ (Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 105)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / LhCtINVERT / ZA
      EXTERNAL CtLhZFXL
      DATA TEM, RER / D1, 1E-3 /
      DATA ZLOW, ZHIGH, IWRN2 / -10.0, 1.00002, 0 /
    7 EPS = TEM * RER
      ZA = Z
      IF (Z .LE. ZHIGH .AND. Z .GT. ZLOW) THEN
          XLA = LOG (XMIN) * 1.5
          XLB = 0.00001
          TEM = CtLhZBRNT (CtLhZFXL, XLA, XLB, EPS, IRT)
      Else
        CALL CtLhWARNR (IWRN2, 'Z out of range in CtLhXFRMZ, X set=0.',
     >              'Z', Z, ZLOW, ZHIGH, 1)
        TEM = 0
      EndIf
      CtLhXFRMZ = EXP(TEM)
      RETURN
      END
      FUNCTION CtLhxLi(n,x)
      implicit NONE
      integer NCUT, i,n,m3
      real*8 CtLhxLi,Out,x,pi2by6,zeta3,c1,c2
      real*8 r,xt,L,xln1m
      parameter (m3=8)
      dimension c1(2:m3),c2(2:m3)
      data NCUT/27/
      data c1/0.75,-0.5833333333333333d0,0.454861111111111d0,
     >        -0.3680555555555555d0,0.3073611111111111d0,
     >        -0.2630555555555555d0,0.2294880243764172d0/
      data c2/-0.5d0,0.5d0,-0.4583333333333333d0,0.416666666666666d0,
     >        -0.3805555555555555d0,0.35d0,-0.3241071428571428d0/
      data zeta3,pi2by6 /1.20205690315959d0,1.64493406684823d0/
      L=0.0
      i=0
      r=1.0
      if (abs(x).gt.r) then
        PRINT *,'Li: x out of range (-1,1) , x=',x
        STOP
      endif
      if (n.lt.0) then
       PRINT *,'Polylogarithm Li undefined for n=',n
       STOP
      elseif (n.eq.0) then
       Out=x/(1d0-x)
      elseif (n.eq.1) then
       Out=-dlog(1-x)
      elseif (n.eq.2) then
                                                !Calculate dilogarithm
                                                !separately for x<0.5 and x>0.5
      if (x.ge.(-0.5).and.x.le.0.5) then
         do while(i.le.NCUT)
       	  i=i+1
          r=r*x
          L=L+r/i/i
         enddo
         Out=L
       elseif (x.eq.0) then
         Out=0d0
       elseif(x.gt.0.5) then !n.eq.2,x>0.5
         xt = 1.0-x
         L = pi2by6 - dlog(x)*dlog(xt)
         do while(i.le.NCUT)
          i=i+1
          r=r*xt
          L=L-r/i/i
         enddo
         Out=L
       elseif (x.lt.(-0.5)) then
         xt=-x/(1d0-x)
         L=-0.5*dlog(1-x)**2
         do while (i.le.NCUT)
          i=i+1
          r=r*xt
          L=L-r/i/i
         enddo
         Out=L
       endif
      elseif (n.eq.3.and.x.ge.0.8) then !use the expansion of Li3 near x=1
       L=zeta3+pi2by6*dlog(x)
       xt=(1d0-x)
       xln1m=dlog(xt)
       do i=2,m3
        L=L+(c1(i)+c2(i)*xln1m)*xt**i
       enddo
       Out=L
      else !n>3 or x=3,x<0.8
         do while(i.le.NCUT)
          i=i+1
          r=r*x
          L=L+r/dble(i)**dble(n)
         enddo
         Out=L
      endif
      CtLhxLi=Out
      End ! CtLhxLi
      FUNCTION CtLhZBRLAM (WLLN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtTRNCOM / VMULM, JRDR, N, N1
      WMULM = EXP (WLLN)
      TEM1 = 1./ CtLhALPQCD(JRDR, N1, WMULM, I)
      TEM2 = 1./ CtLhALPQCD(JRDR, N,  VMULM, I)
      CtLhZBRLAM = TEM1 - TEM2
      END
      FUNCTION CtLhZBRNT(FUNC, X1, X2, TOL, IRT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)
      external func
      IRT = 0
      TOL = ABS(TOL)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        PRINT *, 'Root must be bracketed for CtLhZBRNT. Set = 0'
        IRT = 1
        CtLhZBRNT=0.
        RETURN
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          CtLhZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      PRINT *, 'CtLhZBRNT exceeding maximum iterations.'
      IRT = 2
      CtLhZBRNT=B
      RETURN
      END
      FUNCTION CtLhZFRMX (XX)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 105)
      COMMON / LhCtXXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      DATA IWRN1, HUGE, TINY / 0, 1.E35, 1.E-35 /
      F(X) = (XCR-XMIN) * LOG (X/XMIN) + LOG (XCR/XMIN) * (X-XMIN)
      D(X) = (XCR-XMIN) / X          + LOG (XCR/XMIN)
      X = XX
      IF (X .GE. XMIN) THEN
         TEM = F(X) / F(D1)
      ElseIF (X .GE. D0) THEN
         X = MAX (X, TINY)
         TEM = F(X) / F(D1)
      Else
         CALL CtLhWARNR(IWRN1, 'X out of range in CtLhZFRMX'
     >             , 'X', X, TINY, HUGE, 1)
         TEM = 99.
         STOP
      EndIf
      CtLhZFRMX = TEM
      RETURN
      ENTRY CtLhDZDX (XX)
      X = XX
      IF (X .GE. XMIN) THEN
         TEM = D(X) / F(D1)
      ElseIF (X .GE. D0) THEN
         X = MAX (X, TINY)
         TEM = D(X) / F(D1)
      Else
         CALL CtLhWARNR(IWRN1, 'X out of range in CtLhDZDX '
     >             , 'X', X, TINY, HUGE, 1)
         TEM = 99.
         STOP
      EndIf
      CtLhDZDX = TEM
      RETURN
      END
      FUNCTION CtLhZFXL (XL)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / LhCtINVERT / ZA
      X = EXP(XL)
      TT = CtLhZFRMX (X) - ZA
      CtLhZFXL = TT
      RETURN
      END


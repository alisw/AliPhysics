      SUBROUTINE PHOENE(MPASQR,MCHREN,BETA,IDENT)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays calculation  of photon ENErgy
C.              fraction
C.
C.    Purpose:  Subroutine  returns  photon  energy fraction (in (parent
C.              mass)/2 units) for the decay bremsstrahlung.
C.
C.    Input Parameters:  MPASQR:  Mass of decaying system squared,
C.                       XPHCUT:  Minimum energy fraction of photon,
C.                       XPHMAX:  Maximum energy fraction of photon.
C.
C.    Output Parameter:  MCHREN:  Renormalised mass squared,
C.                       BETA:    Beta factor due to renormalisation,
C.                       XPHOTO:  Photon energy fraction,
C.                       XF:      Correction factor for PHOFAC.
C.
C.    Author(s):  S. Jadach, Z. Was               Created at:  01/01/89
C.                B. van Eijk                     Last Update: 26/03/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION MPASQR,MCHREN,BIGLOG,BETA,DATA
      INTEGER IWT1,IRN,IWT2
      REAL*8 PRSOFT,PRHARD,PHORAN,PHOFAC
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL*8 PNEUTR
      INTEGER IDENT
      REAL*8 PHOCHA
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL*8 ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      REAL*8 PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER IREP
      REAL*8 PROBH,CORWT,XF
      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
      LOGICAL INTERF,ISEC,IFTOP
      REAL*8 FINT,FSEC
      COMMON /PHOKEY/ FSEC,FINT,INTERF,ISEC,IFTOP
C--
      IF (XPHMAX.LE.XPHCUT) THEN
        XPHOTO=0.0D0
        RETURN
      ENDIF
C--   Probabilities for hard and soft bremstrahlung...
      MCHREN=4.D0*MCHSQR/MPASQR/(1.D0+MCHSQR/MPASQR)**2
      BETA=SQRT(1.D0-MCHREN)
      BIGLOG=LOG(MPASQR/MCHSQR*(1.D0+BETA)**2/4.D0*
     &          (1.D0+MCHSQR/MPASQR)**2)
      PRHARD=ALPHA/PI/BETA*BIGLOG*(LOG(XPHMAX/XPHCUT)-.75D0+XPHCUT/
     &XPHMAX-.25D0*XPHCUT**2/XPHMAX**2)
      PRHARD=PRHARD*PHOCHA(IDENT)**2*FINT*FSEC
      IF (IREP.EQ.0) PROBH=0.D0
      PRHARD=PRHARD*PHOFAC(0)
      PROBH=PRHARD
      PRSOFT=1.D0-PRHARD
C--
C--   Check on kinematical bounds
      IF (PRSOFT.LT.0.1D0) THEN
        DATA=PRSOFT
        CALL PHOERR(2,'PHOENE',DATA)
      ENDIF
      IF (PHORAN(IWT1).LT.PRSOFT) THEN
C--
C--   No photon... (ie. photon too soft)
        XPHOTO=0.D0
      ELSE
C--
C--   Hard  photon... (ie.  photon  hard enough).
C--   Calculate  Altarelli-Parisi Kernel
   10   XPHOTO=EXP(PHORAN(IRN)*LOG(XPHCUT/XPHMAX))
        XPHOTO=XPHOTO*XPHMAX
        IF (PHORAN(IWT2).GT.((1.D0+(1.D0-XPHOTO/XPHMAX)**2)/2.D0)) 
     &                            GOTO 10
      ENDIF
C--
C--   Calculate parameter for PHOFAC function
      XF=4.D0*MCHSQR*MPASQR/(MPASQR+MCHSQR-MNESQR)**2
      RETURN
      END

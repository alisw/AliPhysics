C...This file contains a complete listing of all PYTHIA 
C...commonblocks, and additionally some recommended other 
C...declarations. You may copy this to the top of your
C...mina program and then eliminate unnecessary parts.
C  Jun 19 1998 P.Murat(CDF): add implicit for integers
C  Mar 18 1998 P.Murat: remove implicits, add explicit type definition
c                       for all the variables
C-----------------------------------------------------------------
 
C...All real arithmetic in double precision.
c      IMPLICIT DOUBLE  PRECISION(A-H, O-Z)
c      implicit integer (i-n)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers
C...(left- and righthanded SUSY, excited fermions).
      integer  KSUSY1, KSUSY2, KEXCIT
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)

C...Commonblocks.
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer       N,NPAD,K
      real*8        P,V
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer       MSTU, MSTJ
      real*8        PARU, PARJ
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      integer       KCHG
      real*8        PMAS, PARF, VCKM
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer       MDCY, MDME, KFDP
      real*8        BRAT
C...Particle names
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
C...Random number generator information.
      COMMON/PYDATR/MRPY(6),RRPY(100)
      integer       MRPY
      real*8        RRPY
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      integer       MSEL, MSELPD, MSUB, KFIN
      real*8        CKIN
C...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      integer       MSTP, MSTI
      real*8        PARP, PARI
C...Internal variables.
      COMMON/PYINT1/MINT(400),VINT(400)
      integer       mint
      real*8        vint
C...Process information.
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      integer       ISET, KFPR, ICOL
      real*8        COEF
C...Parton distributions and cross sections.
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      integer       ISIG
      real*8        XSFX, SIGH
C...Resonance width and secondary decay treatment.
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      integer       MWID
      real*8        WIDS
C...Generation and cross section statistics.
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      integer       NGENPD, NGEN
      real*8        XSEC
C...Process names.
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
C...Total cross sections.
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      real*8        SIGT
C...Photon parton distributions: total and valence only.
      COMMON/PYINT8/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6), 
     &XPDIR(-6:6) 
      real*8        XPVMD, XPANL, XPANH, XPBEH, XPDIR
      COMMON/PYINT9/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6) 
      real*8        VXPVMD, VXPANL, VXPANH, VXPDGM
C...Setting up user-defined processes.
      COMMON/PYUPPR/NUP,KUP(20,7),NFUP,IFUP(10,2),PUP(20,5),Q2UP(0:10)
      integer       NUP, KUP, NFUP, IFUP
      real*8        PUP, Q2UP
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      integer       IMSS
      real*8        RMSS
C...Supersymmetry mixing matrices.
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      real*8        ZMIX, UMIX, VMIX, SMZ, SMW, SFMIX, ZMIXI,UMIXI,VMIXI
C...Parameters for Gauss integration of supersymmetric widths.
      COMMON/PYINTS/XXM(20)
      real*8        XXM
C...Histogram information.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      integer       IHIST, INDX
      real*8        BIN

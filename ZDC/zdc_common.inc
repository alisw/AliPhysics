*
* Variables for beam divergency and crossing angle (important for ZDC !)
*
      COMMON/BMDIV/BMDIV,BMCRA,IFLCR
      REAL BMDIV,BMCRA
      INTEGER IFLCR
*
*
* Variables relative to size (and offset) of the beam
*
      COMMON /BMTYPE/ BEAMTY,EPERN
      INTEGER BEAMTY
      COMMON /INBEAM/ FX,FY,SIGX,SIGY
*
*
* Flags and arrays for the generation of Fermi momentum of the spectators
*
      COMMON/CFERMI/IFERMI,
     +            PROBINTP(0:200),PROBINTN(0:200),PP(0:200),DDP(3)
*
*
* Geometrical quantities relative to ZDC
*
      COMMON /ZDC_CGEOM/
     + HDZN(3), HDZP(3), FIZN(3), FIZP(3), GRZN(3), GRZP(3),
     + NCEN(3), NCEP(3), ZNPOS(3), ZPPOS(3)
*
*
* Steering flags and variables for HIJING/VENUS
*
      COMMON/EVGN/IHIJ,IHIJF,IHIJSP,IVNU,IVNUF,IVNSP
      COMMON/EVVAR/BIMP,PARTI
*
*
* Files for input/output
*
      CHARACTER*60 FILEA,FILEV,FILEH
      COMMON/FIL/FILEA,FILEV,FILEH
      INTEGER FILEI(15),FILEJ(15),FILEK(15)
      EQUIVALENCE(FILEA,FILEI(1))
      EQUIVALENCE(FILEV,FILEJ(1))
      EQUIVALENCE(FILEH,FILEK(1))
*
*
* Variables for hit storage
*
* Set ZDC
*
      PARAMETER(NVZN=2)
      PARAMETER(NVZP=2)
      COMMON /ZDC_HITS/ ZSET,ZN1,ZP1,CHNMSZN(NVZN),CHNMSZP(NVZP),
     +                  NBITSZN(NVZN),NBITSZP(NVZP),IDTYPZN,IDTYPZP,
     +                  NWHI,NWDI,ITOLD
      CHARACTER*4 ZSET,ZN1,ZP1
      CHARACTER*4 CHNMSZN,CHNMSZP
      INTEGER NBITSZN,NBITSZP
      INTEGER IDTYPZN,IDTYPZP
*
* Detectors ZN1,ZP1
*
      PARAMETER(NHZ=7)
      PARAMETER(NZPTX=2,NZPTY=2,NZNTX=2,NZNTY=2)
      COMMON /ZDC_DET/ CHNAMH(NHZ),NBITSH(NHZ),ORIG(NHZ),FACT(NHZ),
     +                 HITS(NHZ),HITSBUPP(NHZ),HITSBUPN(NHZ),
     +                 NHSUM
      CHARACTER*4 CHNAMH
      INTEGER NBITSH
      REAL ORIG,FACT
      REAL HITS
      REAL HITSBUPP,HITSBUPN
      INTEGER NHSUM
*
* Variables for CWN
*
      COMMON/ZNTVAR/NPA,ENE,SFL,
     +XZP,YZP,
     +EZP(NZPTX,NZPTY),NPP(NZPTX,NZPTY),
     +XZN,YZN,
     +EZN(NZNTX,NZNTY),NPN(NZNTX,NZNTY)
*
      INTEGER NPA,NPP,NPN,SFL
      REAL ENE,XZP,YZP,EZP,XZN,YZN,EZN
*
*
* Flag for showering of primary particles
*
      COMMON/SECFL/ISECF
*
*
* Arrays for Cerenkov photons
*
      PARAMETER (NBETA=4, NALFA=90, NBE=18)
      COMMON /TABL/ TABLE(NBETA,NALFA,NBE)
*

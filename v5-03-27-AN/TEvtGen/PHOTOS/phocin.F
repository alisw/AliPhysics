      SUBROUTINE PHOCIN
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton Common INitialisation
C.
C.    Purpose:  Initialisation of parameters in common blocks.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  Commons /PHOLUN/, /PHOPHO/, /PHOCOP/, /PHPICO/
C.                                and /PHSEED/.
C.
C.    Author(s):  B. van Eijk                     Created at:  26/11/89
C.                Z. Was                          Last Update: 10/08/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      LOGICAL QEDRAD
      COMMON/PHOQED/QEDRAD(NMXHEP)
      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN
      REAL*8 ALPHA,XPHCUT
      COMMON/PHOCOP/ALPHA,XPHCUT
      REAL*8 PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      INTEGER ISEED,I97,J97
      REAL*8 URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
      INTEGER PHOMES
      PARAMETER (PHOMES=10)
      INTEGER STATUS
      COMMON/PHOSTA/STATUS(PHOMES)
      LOGICAL INTERF,ISEC,IFTOP
      REAL*8 FINT,FSEC
      COMMON /PHOKEY/ FSEC,FINT,INTERF,ISEC,IFTOP
      INTEGER INIT,I
      SAVE INIT
      DATA INIT/ 0/
C--
C--   Return if already initialized...
      IF (INIT.NE.0) RETURN
      INIT=1
C--
C--   Preset switch  for  photon emission to 'TRUE' for each particle in
C--   /PH_HEPEVT/, this interface is needed for KORALB and KORALZ...
      DO 10 I=1,NMXHEP
   10 QEDRAD(I)=.TRUE.
C--
C--   Logical output unit for printing of PHOTOS error messages
      PHLUN=6
C--
C--   Set cut parameter for photon radiation
      XPHCUT=0.01D0
C--
C--   Define some constants
      ALPHA=0.00729735039D0
      PI=3.14159265358979324D0
      TWOPI=6.28318530717958648D0
C--
C--   Default seeds Marsaglia and Zaman random number generator
      ISEED(1)=1802
      ISEED(2)=9373
C--
C--   Iitialization for extra options
C--   (1)
C--   Interference weight for two body symmetric channels only. 
      INTERF=.TRUE.
C--   (2)
C--   Second order - double photon switch
      ISEC=.TRUE.
C--   (3)
C--   Emision in the hard process g g (q qbar) --> t tbar 
C--                                 t          --> W b
      IFTOP=.TRUE.
C--
C--   further initialization done automatically
      IF (INTERF) THEN
C--   best choice is if FINT=2**N where N+1 is maximal number 
C--   of charged daughters
C--   see report on overweihted events
        FINT=2.0D0
      ELSE
        FINT=1.0D0
      ENDIF
C--   Initialise status counter for warning messages
      DO 20 I=1,PHOMES
   20 STATUS(I)=0
      RETURN
      END

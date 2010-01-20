      FUNCTION PHORAN(IDUM)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays RANdom number generator based
C.              on Marsaglia Algorithm
C.
C.    Purpose:  Generate  uniformly  distributed  random numbers between
C.              0 and 1.  Super long period:  2**144.  See also:
C.              G. Marsaglia and A. Zaman,  FSU-SCR-87-50,  for seed mo-
C.              difications  to  this version  see:  F. James DD-Report,
C.              November 1988.  The generator  has  to be initialized by
C.              a call to PHORIN.
C.
C.    Input Parameters:   IDUM (integer dummy)
C.
C.    Output Parameters:  Function value
C.
C.    Author(s):  B. van Eijk, G. Marsaglia and   Created at:  27/09/89
C.                A. Zaman                        Last Update: 27/09/89
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PHORAN
      INTEGER IDUM
      INTEGER ISEED,I97,J97
      REAL*8 URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
   10 PHORAN=URAN(I97)-URAN(J97)
      IF (PHORAN.LT.0.D0) PHORAN=PHORAN+1.D0
      URAN(I97)=PHORAN
      I97=I97-1
      IF (I97.EQ.0) I97=97
      J97=J97-1
      IF (J97.EQ.0) J97=97
      CRAN=CRAN-CDRAN
      IF (CRAN.LT.0.D0) CRAN=CRAN+CMRAN
      PHORAN=PHORAN-CRAN
      IF (PHORAN.LT.0.D0) PHORAN=PHORAN+1.D0
      IF (PHORAN.LE.0.D0) GOTO 10
      RETURN
      END

      SUBROUTINE PHORIN
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation  in decays RANdom number generator init
C.
C.    Purpose:  Initialse PHORAN  with  the user  specified seeds in the
C.              array ISEED.  For details  see also:  F. James  CERN DD-
C.              Report November 1988.
C.
C.    Input Parameters:   ISEED(*)
C.
C.    Output Parameters:  URAN, CRAN, CDRAN, CMRAN, I97, J97
C.
C.    Author(s):  B. van Eijk and F. James        Created at:  27/09/89
C.                                                Last Update: 22/02/90
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION DATA
      REAL*8 S,T
      INTEGER I,IS1,IS2,IS3,IS4,IS5,J
      INTEGER ISEED,I97,J97
      REAL*8 URAN,CRAN,CDRAN,CMRAN
      COMMON/PHSEED/ISEED(2),I97,J97,URAN(97),CRAN,CDRAN,CMRAN
C--
C--   Check value range of seeds
      IF ((ISEED(1).LT.0).OR.(ISEED(1).GE.31328)) THEN
        DATA=ISEED(1)
        CALL PHOERR(8,'PHORIN',DATA)
      ENDIF
      IF ((ISEED(2).LT.0).OR.(ISEED(2).GE.30081)) THEN
        DATA=ISEED(2)
        CALL PHOERR(9,'PHORIN',DATA)
      ENDIF
C--
C--   Calculate Marsaglia and Zaman seeds (by F. James)
      IS1=MOD(ISEED(1)/177,177)+2
      IS2=MOD(ISEED(1),177)+2
      IS3=MOD(ISEED(2)/169,178)+1
      IS4=MOD(ISEED(2),169)
      DO 20 I=1,97
        S=0.D0
        T=0.5D0
        DO 10 J=1,24
          IS5=MOD (MOD(IS1*IS2,179)*IS3,179)
          IS1=IS2
          IS2=IS3
          IS3=IS5
          IS4=MOD(53*IS4+1,169)
          IF (MOD(IS4*IS5,64).GE.32) S=S+T
   10   T=0.5D0*T
   20 URAN(I)=S
      CRAN=362436.D0/16777216.D0
      CDRAN=7654321.D0/16777216.D0
      CMRAN=16777213.D0/16777216.D0
      I97=97
      J97=33
      RETURN
      END

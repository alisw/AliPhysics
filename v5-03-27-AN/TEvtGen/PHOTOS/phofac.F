      FUNCTION PHOFAC(MODE)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays control FACtor
C.
C.    Purpose:  This is the control function for the photon spectrum and
C.              final weighting.  It is  called  from PHOENE for genera-
C.              ting the raw photon energy spectrum (MODE=0) and in PHO-
C.              COR to scale the final weight (MODE=1).  The factor con-
C.              sists of 3 terms.  Addition of  the factor FF which mul-
C.              tiplies PHOFAC for MODE=0 and divides PHOFAC for MODE=1,
C.              does not affect  the results for  the MC generation.  An
C.              appropriate choice  for FF can speed up the calculation.
C.              Note that a too small value of FF may cause weight over-
C.              flow in PHOCOR  and will generate a warning, halting the
C.              execution.  PRX  should  be  included for repeated calls
C.              for  the  same event, allowing more particles to radiate
C.              photons.  At  the  first  call IREP=0, for  more  than 1
C.              charged  decay  products, IREP >= 1.  Thus,  PRSOFT  (no
C.              photon radiation  probability  in  the  previous  calls)
C.              appropriately scales the strength of the bremsstrahlung.
C.
C.    Input Parameters:  MODE, PROBH, XF
C.
C.    Output Parameter:  Function value
C.
C.    Author(s):  S. Jadach, Z. Was               Created at:  01/01/89
C.                B. van Eijk                     Last Update: 13/02/90
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 PHOFAC,FF,PRX
      INTEGER MODE
      INTEGER IREP
      REAL*8 PROBH,CORWT,XF
      COMMON/PHOPRO/PROBH,CORWT,XF,IREP
      SAVE PRX,FF
      DATA PRX,FF/ 0.D0, 0.D0/
      IF (MODE.EQ.0) THEN
        IF (IREP.EQ.0) PRX=1.D0
        PRX=PRX/(1.D0-PROBH)
        FF=1.D0
C--
C--   Following options are not considered for the time being...
C--   (1) Good choice, but does not save very much time:
C--       FF=(1.0D0-SQRT(XF)/2.0D0)/(1.0+SQRT(XF)/2.0D0)
C--   (2) Taken from the blue, but works without weight overflows...
C--       FF=(1.D0-XF/(1-(1-SQRT(XF))**2))*(1+(1-SQRT(XF))/SQRT(1-XF))/2
        PHOFAC=FF*PRX
      ELSE
        PHOFAC=1.D0/FF
      ENDIF
      END

      SUBROUTINE GCINIT
      PARAMETER (NWGEANT=3000000)
      COMMON / GCBANK / GMEM(NWGEANT)
*
* Initialise geant system with a default dimension of NWGEANT
*
      CALL GZEBRA(NWGEANT)
      CALL GINIT()
      CALL GZINIT()
*
      END

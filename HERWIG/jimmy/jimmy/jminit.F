*CMZ :          12/05/95  01.10.55  by  Jonathan Butterworth
*-- Author : JMB 03/12/92
      SUBROUTINE JMINIT
C --------------------------------------------------------------------------
C Purpose:  Run initialisation specific to Jimmy. Call from HERWIG
C           run initialisation (HWEINI).
C --------------------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"


      INTEGER I
      EXTERNAL DFUN, DVNOPT
      DOUBLE PRECISION PHADTMP, SMIN

      NLOST=0
      TOTSCAT=0

C     If MI are not switched on, do nothing.
      IF (MSFLAG.EQ.0) RETURN 



C --  Title page.
      WRITE (JMOUT,1001)


C     If MI are not implemented for this process, print warning and set
C     MSFLAG to 0
      IF ((ABS(IPRO).LT.10).OR.(ABS(IPRO).GE.50)) THEN
         
         WRITE(JMOUT,1006) 
     &       'JMINIT: Multiparton interactions make no sense'
     &      ,'        for this process: Turning them off.'
         MSFLAG=0
         RETURN
         
      ELSE IF (ABS(IPRO).EQ.15) THEN
         
c     MI fine.
         WRITE(JMOUT,1006) 'Multiparton interactions being used. '

         IF (JMUEO.EQ.0) THEN
            PTJIM=PTMIN
         ELSE 
            JMUEO=1
            IF (2.0*PTJIM.GT.PTMIN) THEN
               WRITE(6,*)
     $              'JMINIT: WARNING. PTJIM AND PTMIN ARE TOO SIMILAR.'
               WRITE(6,*)
     $              'JMINIT: RECOMMEND YOU USE JMUEO=0'
            ENDIF

            WRITE(JMOUT,1006)
     &      'JMINIT: Multiparton interactions being used in the '
            WRITE(JMOUT,1006)
     &      '        underlying event under high ET jets'
            
         ENDIF

      ELSE

c     MI approx.
         WRITE(JMOUT,1006)
     &      'JMINIT: Multiparton interactions being used in the '
         WRITE(JMOUT,1006)
     &     '         underlying event only.'
         WRITE(JMOUT,1006)
     &     '         Main process is not eikonalised. '

         JMUEO=2

      ENDIF

      WRITE(JMOUT,1008) PTJIM

      WRITE(JMOUT,1009) JMV2,JMU2

C --  Set up the radii appropriately (unless they have been set by the user).
      IF (JMV2.EQ.0) JMV2=JMRAD(IPART1)
      IF (JMU2.EQ.0) JMU2=JMRAD(IPART2)
      IF (JMV2.EQ.0) THEN
        WRITE(JMOUT,8001) 1
	WRITE(JMOUT,*) IPART1,IPART2
	WRITE(JMOUT,*) JMRAD
        JMV2=JMRAD(127)
      ENDIF
      IF (JMU2.EQ.0) THEN
        WRITE(JMOUT,8001) 2
        JMU2=JMRAD(127)
      ENDIF

      WRITE(JMOUT,1009) JMV2,JMU2

C --  Decide whether the energy of 0,1 or 2 of the beams is varying.
C     This assumes that if the beam is an electron, we are dealing with
C     a spectrum of photon energies, but if it is a photon it
C     is monoenergetic. All other beams are monoenergetic hadrons.
      JCMVAR=0
      PHADTMP=PHAD
      PHAD=1.0
      IF (ABS(IDPDG(IPART1)).EQ.11) THEN
         JCMVAR=JCMVAR+1
         PHAD=PHAD*PHADTMP
      ELSE IF (IDPDG(IPART1).EQ.22) THEN
         PHAD=PHAD*PHADTMP
      ENDIF
      IF (ABS(IDPDG(IPART2)).EQ.11) THEN
         JCMVAR=JCMVAR+1
         PHAD=PHAD*PHADTMP
      ELSE IF (IDPDG(IPART2).EQ.22) THEN
         PHAD=PHAD*PHADTMP
      ENDIF

C     Calculate cross section(s)
      WRITE(JMOUT,1005)

c     Minimum allowed Z=S(hadhad)/S(beambeam).
      IF (JCMVAR.EQ.0) THEN
         JMZMIN=1.0
      ELSE
         CALL HWEGAS(SMIN)
         IF (JCMVAR.EQ.1) THEN
            JMZMIN=MAX(YWWMIN,SMIN/
     &           (2.0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2)))
         ELSE IF (JCMVAR.EQ.2) THEN
C     Both particle energies are varying.
            JMZMIN=MAX(YWWMIN**2,SMIN/
     &           (2.0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2)))
         ENDIF
      ENDIF

      CALL JMGRID

c     Restore the photon/hadron probability.
      PHAD=PHADTMP

      IF (JMBUG.GT.0) THEN
        WRITE(JMOUT,1006) 'Z values'
        DO I = 1, NPSIMP
          WRITE(JMOUT,1007) JMARRY(1,I)
        ENDDO
        WRITE(JMOUT,1006) 'UnEikonalised xsec'
        DO I = 1, NPSIMP
          WRITE(JMOUT,1007) JMARRY(2,I)*GEV2NB
        ENDDO
        WRITE(JMOUT,1006) 'Eikonalised xsec'
        DO I = 1, NPSIMP
          WRITE(JMOUT,1007) JMARRY(3,I)*GEV2NB
        ENDDO
      ENDIF

      RETURN

 1001 FORMAT( 8X,'----------------------------------------------'/,
     &        8X,'                J  I  M  M  Y               '/,
     &        8X,'   Inclusion of  Eikonal multiple scattering '/,
     &        8X,'                Version 4.2            '/,
     &        8X,'    Requires HERWIG 6.510 or above   '/,
     &        8X,'----------------------------------------------')
 1005 FORMAT( 8X,'       GENERATING CROSS-SECTION GRID         ')
 1006 FORMAT(8X,A)
 1007 FORMAT(F20.5)
 1008 FORMAT(8X,'MINIMUM PT FOR SECONDARY SCATTERS=',G8.2,'GeV')
 1009 FORMAT(8X,'INVERSE RADII OF COLLIDING PARTICLES(1,2)=',G8.2,' ',g8
     $     .2,'GeV')

 8001 FORMAT(8X,'JIMMY:NO RADIUS FOR PARTICLE ',I2,':USING PHOTON')

 999  CONTINUE
      END

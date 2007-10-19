*$ CREATE CRNKVP.FOR
*COPY CRNKVP
*
*=== Crnkvp ===========================================================*
*
      SUBROUTINE CRNKVP ( MREG, NEWREG, DDNEAR, LEMAGN, DEDXCK )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1997-2003      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     CeReNKoV Production:                                             *
*                                                                      *
*     Created on 21 september 1997 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 30-oct-98     by    Alfredo Ferrari               *
*                                                                      *
*    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    *
*    !!!!!! Important: Etrack, Ptrack, Atrack,  are supposed !!!!!!    *
*    !!!!!! to be the INITIAL values (not yet updated at the !!!!!!    *
*    !!!!!! end of the step).                                !!!!!!    *
*    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    *
*                                                                      *
*    dN/dEdx = alpha^2 z^2 / [ r_e m_e c^2 ] [ 1 - 1 / (beta n)^2 ]    *
*                                                                      *
*    Note the wavelength is always assumed to be the vacuum one        *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(MULBOU)'
      INCLUDE '(OPPHCM)'
      INCLUDE '(OPPHST)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'
*
      PARAMETER ( CSNPRN = 100.D+00 * CSNNRM )
      PARAMETER ( CRNKC0 = ALPFSC * ALPFSC / RCLSEL / AMELCT )
      LOGICAL LEMAGN, LMEQNT, LNWSUB, LSMPAN, LDUMMY, LBXRFL
*  Statement functions:
*     INCLUDE '(DORTSF)'
      INCLUDE '(UNRTSF)'
*
*  No magnetic field end step Newreg check implemented for the moment:
*
*
*      IF ( LEMAGN )
*     &   CALL FLABRT ( 'CRNKVP', ' STOP:LEMAGN-NOT-YET-IMPLEMENTED' )
*  No change of lattice check implemented for the moment:
*      IF ( LT1TRK .NE. LT2TRK ) CALL FLABRT ( 'CRNKVP',
*     &    ' STOP:LT1TRK.NE.LT2TRK-NOT-YET-IMPLEMENTED' )
      DEDXCK = ZERZER
      MMAT   = MEDIUM (MREG)
      IF ( JTRACK .GE. -6 ) THEN
         AMPRTC = AM (JTRACK)
      ELSE
         AMPRTC = AMNHEA (-JTRACK)
      END IF
*  +-------------------------------------------------------------------*
*  |  Find the maximum refractive index over the interval allowed
*  |  for Cerenkov production:
      IF ( RMXCER (MMAT) .LT. ZERZER ) THEN
         RMXCER (MMAT) = RFNDMX (MMAT)
      END IF
*  |
*  +-------------------------------------------------------------------*
      BETNCR = PTRACK / ETRACK * RMXCER (MMAT)
*  Maximum beta x n already below 1: no chance to emit any photon
      IF ( BETNCR .LE. ONEONE ) RETURN
      LMEQNT = MTRACK .EQ. NTRACK
      DTRCKT = ZERZER
*  +-------------------------------------------------------------------*
*  |  Loop on energy deposition sub-steps:
      DO 300 ID = 1, MTRACK
         DTRCKT = DTRCKT + DTRACK (ID)
  300 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      TTRCKT = ZERZER
*  +-------------------------------------------------------------------*
*  |  Loop on track sub-steps (be careful with magnetic field they are
*  |  not the straight distance between the end points):
      DO 500 IT = 1, NTRACK
         TTRCKT = TTRCKT + TTRACK (IT)
  500 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      CRVCRR = CTRACK / TTRCKT
      IF ( JTRACK .GE. -6 ) THEN
         CRNKCS = CRNKC0 * ( DBLE ( ICHRGE (JTRACK) ) )**2 * OPSNMX
      ELSE
         CRNKCS = CRNKC0 * ( DBLE ( ICHEAV(-JTRACK) ) )**2 * OPSNMX
      END IF
      BTNFCR = ( BETNCR - ONEONE ) * ( BETNCR + ONEONE ) / BETNCR**2
*  Macroscopic sigma for production (cm^-1)
      SIGMCK = CRNKCS * ( EMXCER (MMAT) - EMNCER (MMAT) ) * BTNFCR
*  Average number of emitted photons:
      NEMPHO = NINT ( CTRACK * SIGMCK )
*  Take a three sigma clearance:
      NEMPHO = NEMPHO + 10 + 3 * NINT ( SQRT ( CTRACK * SIGMCK ) )
      NEMPHO = MIN ( NEMPHO, MOSTCK )
*D     IF ( LSTOPP .GT. MOSTCK - NEMPHO ) WRITE (77,*)
*D    &   ' ###CRNKVP:LSTOPP:',LSTOPP,NEMPHO
*  +-------------------------------------------------------------------*
*  |  Empty the stack if too full: compute and save the normal in case
*  |  it is required
      IF ( LSTOPP .GT. MOSTCK - NEMPHO ) THEN
         LDUMMY = LBXRFL ( NEWREG, MREG, ZERZER, ZERZER, ONEONE )
         CALL KASOPH ( .TRUE. )
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Total wiggled and curved track used up to now:
      TRUSED = ZERZER
      ETINTR = ETRACK
      PTINTR = PTRACK
      TRINTR = ZERZER
      NTRKLD = 0
*  Current lattice number and index inside Mulbou:
      IRLTCR = 0
      LATCCR = MULTTC (IRLTCR)
*  Variables in the current Ntrckr bin:
      NTRKCR = 1
      PTRKCR = PTRACK
      ETRKCR = ETRACK
      ATRKCR = ATRACK
      TRNRSD = TTRACK (NTRKCR) * CRVCRR
*  +-------------------------------------------------------------------*
*  |  Ntrack = Mtrack
      IF ( LMEQNT ) THEN
         DEDXCR = DTRACK (NTRKCR) / TRNRSD
*  |
*  +-------------------------------------------------------------------*
*  |  Ntrack >< Mtrack
      ELSE
         DEDXCR = DTRCKT / CTRACK
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |  Stacking loop for Cerenkov photons:
      NPROD = 0
 1000 CONTINUE
*D        IF ( SIGMCK .LT. ZERZER ) WRITE (77,*)
*D    &      ' ^^^CRNKVP:SIGMCK,BTNFCR',SIGMCK,BTNFCR
*D        IF ( DEDXCR .LT. AZRZRZ ) WRITE (77,*)' ###CRNKVP:',
*D    &  'MMAT,JTRACK,DEDXCR,NTRKCR,DTRACK(NTRKCR),CTRACK,TRNRSD,DTRCKT',
*D    &   MMAT,JTRACK,DEDXCR,NTRKCR,DTRACK(NTRKCR),CTRACK,TRNRSD,DTRCKT
*  |  Empty the stack if too full:
*D        IF ( LSTOPP .EQ. MOSTCK ) WRITE (77,*)
*D    &      ' ###CRNKVP:LSTOPP:',LSTOPP
*  |  +----------------------------------------------------------------*
*  |  |  Empty the stack if too full: compute and save the normal in
*  |  |  case it is required
         IF ( LSTOPP .EQ. MOSTCK ) THEN
            LDUMMY = LBXRFL ( NEWREG, MREG, ZERZER, ZERZER, ONEONE )
            CALL KASOPH ( .TRUE. )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         DSTPRD = - LOG ( ONEONE - FLRNDM (DSTPRD) ) / SIGMCK
         TRUSED = TRUSED + DSTPRD
*dbgD        WRITE (77,*) ' CRNKVP:DSTPRD,TRUSED,CTRACK,TRNRSD',
*dbgD    &                         DSTPRD,TRUSED,CTRACK,TRNRSD
         IF ( TRUSED .GE. CTRACK ) GO TO 7000
*  |  +----------------------------------------------------------------*
*  |  |  Sub-step loop:
 1500    CONTINUE
*D           IF ( NTRKCR .GT. NTRACK ) WRITE (77,*)
*D    &         ' ^^^CRNKVP:NTRKCR,NTRACK',NTRKCR,NTRACK
*  |  |  +-------------------------------------------------------------*
*  |  |  |  The production point is inside the current sub-step:
            IF ( DSTPRD .LE. TRNRSD ) THEN
*  |  |  |  Update the distance from the previous interaction:
               TRINTR = TRINTR + DSTPRD
*  |  |  |  Update the residual distance available in the current
*  |  |  |  sub-step
               TRNRSD = TRNRSD - DSTPRD
               ETRKCR = ETRKCR - DSTPRD * DEDXCR
               LNWSUB = .FALSE.
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |  The production point is outside the current sub-step:
            ELSE
*  |  |  |  Update the distance from the previous interaction:
               TRINTR = TRINTR + TRNRSD
               ETRKCR = ETRKCR - TRNRSD * DEDXCR
               DSTPRD = DSTPRD - TRNRSD
*  |  |  |  Change sub-step:
               NTRKCR = NTRKCR + 1
*  |  |  |  Update the residual distance available in the current
*  |  |  |  sub-step
               TRNRSD = TTRACK (NTRKCR) * CRVCRR
*  |  |  |  Update the current dE/dx:
               IF ( LMEQNT ) DEDXCR = DTRACK (NTRKCR) / TRNRSD
               LNWSUB = .TRUE.
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         IF ( LNWSUB ) GO TO 1500
*  |  |  end sub-step loop:
*  |  +----------------------------------------------------------------*
         PTRKCR = SQRT ( ( ETRKCR - AMPRTC ) * ( ETRKCR + AMPRTC ) )
*  |  +----------------------------------------------------------------*
*  |  |  Update the particle "age" (beta=p/E, dp/dE=E/p):
*  |  |           Etrkcr
*  |  |         /
*  |  |  Dt = - | dE dx/dE E/(pc) = 1/c dx/dE [ Ptintr - Ptrkcr ]
*  |  |         /
*  |  |          Etintr
         IF ( ETINTR - ETRKCR .GT. CSNNRM * ETRKCR ) THEN
*  |  |  Average dE/dx from the previous (fictitious) interaction:
            DEDXAV = ( ETINTR - ETRKCR ) / TRINTR
            ATRKCR = ATRKCR + ( PTINTR - PTRKCR ) / DEDXAV / CLIGHT
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            ATRKCR = ATRKCR + TRINTR * ETRKCR / PTRKCR / CLIGHT
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         PTINTR = PTRKCR
         ETINTR = ETRKCR
         TRINTR = ZERZER
         BETNLD = BETNCR
         BETNCR = PTRKCR / ETRKCR * RMXCER (MMAT)
*dbgD        WRITE (77,*) ' CRNKVP:BETNCR,RMXCER,LSTOPP',
*dbgD    &                         BETNCR,RMXCER(MMAT),LSTOPP
*  |  No longer any chance to emit Cerenkov photons:
         IF ( BETNCR .LE. ONEONE ) GO TO 7000
         BTNFLD = BTNFCR
         BTNFCR = ( BETNCR - ONEONE ) * ( BETNCR + ONEONE ) / BETNCR**2
         FREJE  = BTNFCR / BTNFLD
*D        IF ( FREJE .GT. ONEPLS ) WRITE (77,*) ' ^^^CRNKVP:',
*D    &       'FREJE,BETNCR,BETNLD,DEDXCR,PTRKCR,ETRKCR,PTRACK,ETRACK',
*D    &        FREJE,BETNCR,BETNLD,DEDXCR,PTRKCR,ETRKCR,PTRACK,ETRACK
*  |  Update the macroscopic sigma: beta can only decrease
         SIGMCK = SIGMCK * FREJE
         RNDREJ = FLRNDM (RNDREJ)
         IF ( RNDREJ .GE. FREJE ) GO TO 1000
*  |--<--<--<--<--< Interaction rejected because of the decrease in beta
         RNDREJ = RNDREJ / FREJE
*  |  Sample the emitted photon energy:
         EPHSMP = RNDREJ * ( EMXCER (MMAT) - EMNCER (MMAT) )
     &          + EMNCER (MMAT)
*  |  2pi x freq.:
         OMGSMP = GEVOMG * EPHSMP
*  |  Wavelength:
         WVLSMP = TWOPIP * CLIGHT / OMGSMP
*  |  Compute the refraction index at the right wavelength/frequency:
         RFISMP = FOPTPR ( 1, WVLSMP, OMGSMP, MMAT )
         BETNSM = PTRKCR / ETRKCR * RFISMP
*dbgD        WRITE (77,*) ' CRNKVP:BETNSM,RFISMP,EPHSMP,WVLSMP,OMGSMP',
*dbgD    &                         BETNSM,RFISMP,EPHSMP,WVLSMP,OMGSMP
         IF ( BETNSM .LE. ONEONE ) GO TO 1000
*  |--<--<--<--<--< No chance to emit Cerenkov photons at this frequency
         BTNFSM = ( BETNSM - ONEONE ) * ( BETNSM + ONEONE ) / BETNSM**2
*  |  Compute the quantum efficiency for this emitted energy:
         OPSENS = FOPTSN ( WVLSMP, OMGSMP )
         FREJE  = BTNFSM / BTNFCR * OPSENS / OPSNMX
*D        IF ( FREJE .GT. ONEPLS ) WRITE (77,*)
*D    &      ' ^^^CRNKVP:FREJE,RFISMP,RMXCER(MMAT),OPSENS,MMAT',
*D    &                  FREJE,RFISMP,RMXCER(MMAT),OPSENS,MMAT
         RNDREJ = FLRNDM (RNDREJ)
*dbgD        WRITE (77,*) ' CRNKVP:BTNFSM,BTNFCR,RNDREJ',
*dbgD    &                         BTNFSM,BTNFCR,RNDREJ
         IF ( RNDREJ .GE. FREJE ) GO TO 1000
*  |--<--<--<--<--< Interaction rejected because of the decrease in n
*  |  Eventually a Cerenkov photon is going to be emitted:
*  |  +----------------------------------------------------------------*
*  |  |  New production sub-step:
         IF ( NTRKCR .NE. NTRKLD ) THEN
            UTRKCR = XTRACK (NTRKCR) - XTRACK (NTRKCR-1)
            VTRKCR = YTRACK (NTRKCR) - YTRACK (NTRKCR-1)
            WTRKCR = ZTRACK (NTRKCR) - ZTRACK (NTRKCR-1)
            TUVWCR = SQRT ( UTRKCR**2 + VTRKCR**2 + WTRKCR**2 )
            UTRKCR = UTRKCR / TUVWCR
            VTRKCR = VTRKCR / TUVWCR
            WTRKCR = WTRKCR / TUVWCR
            SINT02 = UTRKCR**2 + VTRKCR**2
            LSMPAN = SINT02 .LT. ANGLSQ
            COSTH0 = WTRKCR
            IF ( .NOT. LSMPAN ) THEN
               SINTH0 = SQRT (SINT02)
               COSPH0 = UTRKCR / SINTH0
               SINPH0 = VTRKCR / SINTH0
            END IF
            NTRKLD = NTRKCR
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         SRNRSD = TUVWCR * TRNRSD / CRVCRR / TTRACK (NTRKCR)
*D        IF ( SRNRSD .LT. ZERZER .OR. SRNRSD .GT. TUVWCR )
*D    &      WRITE (77,*)' ^^^CRNKVP:SRNRSD,TUVWCR,NTRKCR',
*D    &                              SRNRSD,TUVWCR,NTRKCR
         XTRKCR = XTRACK (NTRKCR) - UTRKCR * SRNRSD
         YTRKCR = YTRACK (NTRKCR) - VTRKCR * SRNRSD
         ZTRKCR = ZTRACK (NTRKCR) - WTRKCR * SRNRSD
*  |  +----------------------------------------------------------------*
*  |  |  Put here the check on lattice change:
         IF ( NULTTC .GT. 0 ) THEN
*  |  |  Current (possibly curved in mag. field) non wiggled path:
            TRLATT = TRUSED / CRVCRR
            IF ( TRLATT .LE. TSLTTC (IRLTCR) ) THEN
               LATCCR = MULTTC (IRLTCR)
            ELSE
               DO 2000 IL = IRLTCR + 1, NULTTC
                  IF ( TRLATT .LE. TSLTTC (IL) ) THEN
                     IRLTCR = IL
                     LATCCR = MULTTC (IRLTCR)
                     GO TO 2010
                  END IF
 2000          CONTINUE
*D              CALL FLABRT ( 'CRNKVP', 'STOP:CRNKVP-NO-VALID-LATTICE' )
 2010          CONTINUE
            END IF
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |  Put here the check on the end point in magnetic field:
         IF ( LEMAGN .AND. NTRKCR .EQ. NTRACK .AND. MREG .NE. NEWREG )
     &      THEN
            NEWLSV = NEWLAT
            MLATSV = MLATTC
            NROLD  = NEWREG
            CALL GEOMAG ( XTRKCR, YTRKCR, ZTRKCR, UTRKCR, VTRKCR,
     &                    WTRKCR, NWREG , NROLD , IDSCCK )
            NEWLAT = NEWLSV
            MLATTC = MLATSV
            IF ( NWREG .EQ. NROLD ) GO TO 7000
*  |  |-->-->--> already beyond the real crossing point:
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         LSTOPP = LSTOPP + 1
         NUMOPH = NUMOPH + 1
         IF ( NUMOPH .GT. 1000000000 ) THEN
            MUMOPH = MUMOPH + 1
            NUMOPH = NUMOPH + 1000000000
         END IF
         WOPTPH = WOPTPH + WTRACK
         POPTPH (LSTOPP) = EPHSMP
*  |  +----------------------------------------------------------------*
*  |  |  No easy safe way to reconstruct Dnear along a magnetic field
*  |  |  step:
         IF ( LEMAGN ) THEN
            DONEAR (LSTOPP) = ZERZER
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |  Questionable: ddnear actually refers to the end track point
         ELSE
            DONEAR (LSTOPP) = MAX ( DDNEAR - CTRACK + TRUSED, ZERZER )
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         XOPTPH (LSTOPP) = XTRKCR
         YOPTPH (LSTOPP) = YTRKCR
         ZOPTPH (LSTOPP) = ZTRKCR
         NREGOP (LSTOPP) = MREG
         NLATOP (LSTOPP) = LATCCR
         COSTHE = ONEONE / BETNSM
         SINTHE = SQRT ( ( ONEONE - COSTHE ) * ( ONEONE + COSTHE ) )
         CALL SFECFE ( SINPHI, COSPHI )
         UPRIME = COSPHI * SINTHE
         VPRIME = SINPHI * SINTHE
         WPRIME = COSTHE
         IF ( LSMPAN ) THEN
            TXOPPH (LSTOPP) = UPRIME
            TYOPPH (LSTOPP) = VPRIME
            TZOPPH (LSTOPP) = WPRIME * COSTH0
         ELSE
            TXOPPH (LSTOPP) = UNDOXR ( UPRIME, VPRIME, WPRIME, SINPH0,
     &                                 COSPH0, SINTH0, COSTH0 )
            TYOPPH (LSTOPP) = UNDOYR ( UPRIME, VPRIME, WPRIME, SINPH0,
     &                                 COSPH0, SINTH0, COSTH0 )
            TZOPPH (LSTOPP) = UNDOZR ( UPRIME, VPRIME, WPRIME, SINPH0,
     &                                 COSPH0, SINTH0, COSTH0 )
         END IF
*  |  Set-up the polarization vector: the following for linear pola-
*  |  rization lying in the plane defined by the photon and the
*  |  original particle direction:
         TXPOPP (LSTOPP) = UTRKCR / SINTHE - COSTHE / SINTHE
     &                   * TXOPPH (LSTOPP)
         TYPOPP (LSTOPP) = VTRKCR / SINTHE - COSTHE / SINTHE
     &                   * TYOPPH (LSTOPP)
         TZPOPP (LSTOPP) = WTRKCR / SINTHE - COSTHE / SINTHE
     &                   * TZOPPH (LSTOPP)
         SCADOT = TXOPPH (LSTOPP) * TXPOPP (LSTOPP)
     &          + TYOPPH (LSTOPP) * TYPOPP (LSTOPP)
     &          + TZOPPH (LSTOPP) * TZPOPP (LSTOPP)
         IF ( SCADOT .GT. CSNPRN ) THEN
            WRITE (LUNERR,'(A,1PG23.15,A)')
     &           ' *** Crnkvp: bad polarization',
     &              SCADOT, ' ***'
*D           WRITE (77,'(A,1PG23.15)')
*D    &            ' ^^^Crnkvp: bad polarization',
*D    &              SCADOT
         END IF
         WTOPPH (LSTOPP) = WTRACK
         AGOPPH (LSTOPP) = ATRKCR
         CMPOPP (LSTOPP) = ZERZER
         LOOPPH (LSTOPP) = LTRACK + 1
*
*
*        Hook to TFluka
*
         PXCR =  EPHSMP * TXOPPH (LSTOPP)
         PYCR =  EPHSMP * TYOPPH (LSTOPP)
         PZCR =  EPHSMP * TZOPPH (LSTOPP)
         POX  = TXPOPP(LSTOPP)
         POY  = TYPOPP(LSTOPP)
         POZ  = TZPOPP(LSTOPP)
         CALL pshckp(PXCR, PYCR, PZCR, EPHSMP, XTRKCR, 
     &        YTRKCR , ZTRKCR, ATRKCR, POX, POY, POZ, WTRACK, ITFL)
         NPROD = NPROD + 1
         CALL ustckv(NPROD, MREG, XTRKCR, YTRKCR, ZTRKCR)
*
*
*
*  |  !!!!!! Here Stuprf should be used !!!!!!
         LOUOPP (LSTOPP) = ITFL
         DO 2100 ISPR = 1, MKBMX1
            SPAROK (ISPR,LSTOPP) = SPAUSR (ISPR)
 2100    CONTINUE
         DO 2200 ISPR = 1, MKBMX2
            ISPORK (ISPR,LSTOPP) = ISPUSR (ISPR)
 2200    CONTINUE
*  |  Save the parent features:
         TPROPP (LSTOPP) = ETRACK - AMPRTC
         APROPP (LSTOPP) = ATRKCR
         LPROPP (LSTOPP) = LTRACK
         IPROPP (LSTOPP) = JTRACK
         NPROPP (LSTOPP) = 0
*  |  Update the "current" particle values:
         ETRKCR = ETRKCR - EPHSMP
         DEDXCK = DEDXCK + EPHSMP
         PTRKCR = SQRT ( ( ETRKCR - AMPRTC ) * ( ETRKCR + AMPRTC ) )
         PTINTR = PTRKCR
         ETINTR = ETRKCR
         TRINTR = ZERZER
         BETNLD = BETNCR
         BETNCR = PTRKCR / ETRKCR * RMXCER (MMAT)
*  |  No longer any chance to emit Cerenkov photons:
         IF ( BETNCR .LE. ONEONE ) GO TO 7000
         BTNFLD = BTNFCR
         BTNFCR = ( BETNCR - ONEONE ) * ( BETNCR + ONEONE ) / BETNCR**2
*  |  Update the macroscopic sigma: beta can only decrease
         SIGMCK = SIGMCK * BTNFCR / BTNFLD
      GO TO 1000
*  |
*  +-------------------------------------------------------------------*
 7000 CONTINUE
      RETURN
*=== End of subroutine Crnkvp =========================================*
      END


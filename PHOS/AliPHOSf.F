*CMZ :          15/02/99  14.40.30  by  Federico Carminati
*CMZ :  2.03/01 20/08/98  09.43.00  by  Federico Carminati
*CMZ :  2.03/00 09/07/98  19.21.29  by  Federico Carminati
*CMZ :  2.00/05 25/05/98  14.39.01  by  Federico Carminati
*CMZ :  1.05/06 26/10/95  17.04.53  by  Nick van Eijndhoven (RUU/CERN)
*CMZ :  1.05/00 11/11/94  15.21.30  by  Nick van Eijndhoven (RUU/CERN)
*-- Author :    Nick van Eijndhoven (CERN)   24/09/90
      SUBROUTINE PHOS_INIT
C
C *** DEFINITION OF THE GEOMETRY OF THE PHOS ***
C *** NVE 24-SEP-1990 CERN GENEVA ***
C
C CALLED BY : SXGEOM
C ORIGIN    : NICK VAN EIJNDHOVEN
C
*KEEP,GCKINE.
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      INTEGER       IKINE,ITRA,ISTAK,IVERT,IPART,ITRTYP,NAPART,IPAOLD
      REAL          PKINE,AMASS,CHARGE,TLIFE,VERT,PVERT
C
*KEEP,LUDAT1.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      INTEGER MSTU,MSTJ
      REAL    PARU,PARJ
      SAVE /LUDAT1/
*
*KEEP,SHRUNP.
      COMMON /SHRUNP/ VMAJ,IMIN,NRUN,NEVTOT
*KEEP,SHPHYP.
      COMMON /SHPHYP/ JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,
     +                  JKAC,JKA0,JRHO,JOME,JPHI,JPSI,JDRY
*KEEP,SCPHOS.
        INTEGER MAXCRAD
        PARAMETER (MAXCRAD=100)
        INTEGER PHOSsize,PHOS_Ndiv_magic
        REAL    PHOSflags,PHOScell,PHOSradius,PHOSCPV,
     +          PHOScradlesA,PHOSTXW,PHOSAIR,PHOSFTI,
     +          PHOSextra, PHOSangle
        COMMON /PHOS_PARS/ PHOSflags(9),
     +                           PHOScell(9),PHOSradius,PHOSCPV(9),
     +                           PHOSsize(3), PHOScradlesA,
     +                           PHOSTXW(3),PHOSAIR(3),PHOSFTI(4),
     +                           PHOSextra(9), PHOSangle(MAXCRAD),
     +                           PHOS_Ndiv_magic
*KEND.
C
      IF (PHOSsize(3).GT.MAXCRAD) THEN
         STOP 'PHOS_INIT Increase MAXCRAD, too many cradles!'
      ENDIF

*AZ { SHAKER initialization
      IF( IKINE.EQ.700 ) THEN
         JWEI        = 0
         NDNDY       = NINT(PKINE(4))
         YLIM        = PKINE(5)
         PTLIM       = PKINE(6)
         JWEAK       = 0
         JPI0        = 1
         JETA        = 1
         JPIC        = NINT(PKINE(7))
         JPRO        = NINT(PKINE(7))
         JKAC        = NINT(PKINE(7))
         JKA0        = NINT(PKINE(7))
         JRHO        = NINT(PKINE(7))
         JOME        = NINT(PKINE(7))
         JPHI        = NINT(PKINE(7))
         JPSI        = NINT(PKINE(7))
         JDRY        = NINT(PKINE(7))
         NEVTOT      = 999999999

         MSTU(11)    = 6        ! Shaker output to screen (Fortran channel 6).
         CALL SHINIT            ! Shaker initialization.

         print *,'************************************************'
         print *,'Print some SHAKER parameters.'
         print *,'NDNDY=',NDNDY
         print *,'YLIM=',YLIM
         print *,'PTLIM=',PTLIM
         print *,'JPIC=',JPIC
         print *,'************************************************'
*AZ }

      ENDIF
*
      END
*CMZ :          15/02/99  14.40.30  by  Federico Carminati
*CMZ :  2.03/01 28/07/98  16.12.22  by  Federico Carminati
*-- Author :    Federico Carminati   17/07/98
      SUBROUTINE PHOS_KINE(NT)
*KEEP,GCKINE.
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      INTEGER       IKINE,ITRA,ISTAK,IVERT,IPART,ITRTYP,NAPART,IPAOLD
      REAL          PKINE,AMASS,CHARGE,TLIFE,VERT,PVERT
C
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,BLUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/

*KEEP,SCPHOS.
        INTEGER MAXCRAD
        PARAMETER (MAXCRAD=100)
        INTEGER PHOSsize,PHOS_Ndiv_magic
        REAL    PHOSflags,PHOScell,PHOSradius,PHOSCPV,
     +          PHOScradlesA,PHOSTXW,PHOSAIR,PHOSFTI,
     +          PHOSextra, PHOSangle
        COMMON /PHOS_PARS/ PHOSflags(9),
     +                           PHOScell(9),PHOSradius,PHOSCPV(9),
     +                           PHOSsize(3), PHOScradlesA,
     +                           PHOSTXW(3),PHOSAIR(3),PHOSFTI(4),
     +                           PHOSextra(9), PHOSangle(MAXCRAD),
     +                           PHOS_Ndiv_magic
*KEND.
*
      DIMENSION V_TMP(3),P_TMP(3)
*
      IF(IKINE.EQ.700) THEN
* Shaker genrator          SHOULD BE CHECKED!!!!
         IEVT = IEVT + 1        ! Increament SHAKER event number.
         CALL SHEVNT


***********************************
* Begin:   create list of SHAKER particles (if required)
* ------
* PHOSflags:    YES: X<>0   NO: X=0
* PHOSflags(1) : -----X  Create branch for TObjArray of AliPHOSCradle
*                ----X-  Create file (ftn03 on HP-UX) with list of SHAKER particles (7Mb/event)
         tmp = PHOSflags(1)/(10**2)
         i   = tmp
         i   = (tmp-i)*10

         IF( i.NE.0. ) CALL SHLIST   ! List of SHAKER particles. 7 Mbytes/event

* End of creation
************************************

         DO I=1,N
                                !? IF( K_LUJETS(I,5).EQ.0 )      MPRIMA=MPRIMA+1
            IPART = IPART_FROM_LUJET_TO_GEANT(K(I,2))

            IF( K(I,5).EQ.0 ) THEN
              IFLAG_DONE = 1
            ELSE
              IFLAG_DONE = 0
            ENDIF
            DO JJ=1,3
               V_TMP(JJ)=V(I,JJ)
               P_TMP(JJ)=P(I,JJ)
            ENDDO
            CALL RXSTRAK(IFLAG_DONE,K(I,3),IPART,P_TMP
     $           ,V_TMP,V(I,4),'Primary',NT)
            IF( I.NE.NT ) THEN
               STOP 'GUKINE: Bad thing...'
            ENDIF

            CALL RXKEEP(I)
         ENDDO

      ENDIF
      GOTO 999
*
  999 END
      INTEGER FUNCTION IPART_FROM_LUJET_TO_GEANT(N)
        ! Return GEANT particle number from LUJET particle number.
        ! This code is from file shake005.f :  SUBROUTINE SHTOGL

        IMPLICIT NONE
        INTEGER N

	IF      (N.EQ.0) THEN
          IPART_FROM_LUJET_TO_GEANT = 0
	ELSE IF (N.EQ.22) THEN
	  IPART_FROM_LUJET_TO_GEANT = 1         ! gamma
	ELSE IF (N.EQ.-11) THEN
	  IPART_FROM_LUJET_TO_GEANT = 2
	ELSE IF (N.EQ.11) THEN
	  IPART_FROM_LUJET_TO_GEANT = 3
	ELSE IF (ABS(N).EQ.12) THEN
          IPART_FROM_LUJET_TO_GEANT = 4
	ELSE IF (ABS(N).EQ.14) THEN
	  IPART_FROM_LUJET_TO_GEANT = 4
	ELSE IF (ABS(N).EQ.16) THEN
	  IPART_FROM_LUJET_TO_GEANT = 4
	ELSE IF (N.EQ.-13) THEN
	  IPART_FROM_LUJET_TO_GEANT = 5
	ELSE IF (N.EQ.13) THEN
	  IPART_FROM_LUJET_TO_GEANT = 6
	ELSE IF (N.EQ.111) THEN
	  IPART_FROM_LUJET_TO_GEANT = 7
	ELSE IF (N.EQ.211) THEN
	  IPART_FROM_LUJET_TO_GEANT = 8
	ELSE IF (N.EQ.-211) THEN
	  IPART_FROM_LUJET_TO_GEANT = 9
	ELSE IF (N.EQ.130) THEN
	  IPART_FROM_LUJET_TO_GEANT = 10
	ELSE IF (N.EQ.321) THEN
	  IPART_FROM_LUJET_TO_GEANT = 11
	ELSE IF (N.EQ.-321) THEN
	  IPART_FROM_LUJET_TO_GEANT = 12
	ELSE IF (N.EQ.2112) THEN
	  IPART_FROM_LUJET_TO_GEANT = 13
	ELSE IF (N.EQ.2212) THEN
	  IPART_FROM_LUJET_TO_GEANT = 14
	ELSE IF (N.EQ.-2212) THEN
	  IPART_FROM_LUJET_TO_GEANT = 15
	ELSE IF (N.EQ.310) THEN
	  IPART_FROM_LUJET_TO_GEANT = 16
	ELSE IF (N.EQ.221) THEN
	  IPART_FROM_LUJET_TO_GEANT = 17        ! eta
	ELSE IF (N.EQ.311) THEN                 ! This is (K0) and we set it to
	  IPART_FROM_LUJET_TO_GEANT = 10        ! GEANT K0 short
	ELSE IF (N.EQ.-311) THEN                ! This is (K~0) and we set it to
	  IPART_FROM_LUJET_TO_GEANT = 10        ! GEANT K0 short
	ELSE
          WRITE(*,*) 'Unknown LUJET particle ', N
          stop
	ENDIF
      END

********************************************************************************

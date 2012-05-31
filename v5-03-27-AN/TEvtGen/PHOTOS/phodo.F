      SUBROUTINE PHODO(IP,NCHARB,NEUDAU)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in  decays DOing of KINematics
C.
C.    Purpose:  Starting  from   the  charged  particle energy/momentum,
C.              PNEUTR, photon  energy  fraction and photon  angle  with
C.              respect  to  the axis formed by charged particle energy/
C.              momentum  vector  and PNEUTR, scale the energy/momentum,
C.              keeping the original direction of the neutral system  in
C.              the lab. frame untouched.
C.
C.    Input Parameters:   IP:      Pointer  to   decaying  particle   in
C.                                 /PHOEVT/  and   the   common   itself
C.                        NCHARB:  pointer to the charged radiating
C.                                 daughter in /PHOEVT/.
C.                        NEUDAU:  pointer to the first neutral daughter
C.    Output Parameters:  Common /PHOEVT/, with photon added.
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 27/05/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION PHOAN1,PHOAN2,ANGLE,FI1,FI3,FI4,FI5,TH1,TH3,TH4
      DOUBLE PRECISION PARNE,QNEW,QOLD,DATA
      INTEGER IP,FI3DUM,I,J,NEUDAU,FIRST,LAST
      INTEGER NCHARB
      REAL*8 EPHOTO,PMAVIR,PHOTRI
      REAL*8 GNEUT,PHORAN,CCOSTH,SSINTH,PVEC(4)
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
      DOUBLE PRECISION MCHSQR,MNESQR
      REAL*8 PNEUTR
      COMMON/PHOMOM/MCHSQR,MNESQR,PNEUTR(5)
      DOUBLE PRECISION COSTHG,SINTHG
      REAL*8 XPHMAX,XPHOTO
      COMMON/PHOPHS/XPHMAX,XPHOTO,COSTHG,SINTHG
      REAL*8 PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
C--
      EPHOTO=XPHOTO*PPHO(5,IP)/2.D0
      PMAVIR=SQRT(PPHO(5,IP)*(PPHO(5,IP)-2.D0*EPHOTO))
C--
C--   Reconstruct  kinematics  of  charged particle  and  neutral system
      FI1=PHOAN1(PNEUTR(1),PNEUTR(2))
C--
C--   Choose axis along  z of  PNEUTR, calculate  angle  between x and y
C--   components  and z  and x-y plane and  perform Lorentz transform...
      TH1=PHOAN2(PNEUTR(3),SQRT(PNEUTR(1)**2+PNEUTR(2)**2))
      CALL PHORO3(-FI1,PNEUTR(1))
      CALL PHORO2(-TH1,PNEUTR(1))
C--
C--   Take  away  photon energy from charged particle and PNEUTR !  Thus
C--   the onshell charged particle  decays into virtual charged particle
C--   and photon.  The virtual charged  particle mass becomes:
C--   SQRT(PPHO(5,IP)*(PPHO(5,IP)-2*EPHOTO)).  Construct  new PNEUTR mo-
C--   mentum in the rest frame of the parent:
C--   1) Scaling parameters...
      QNEW=PHOTRI(PMAVIR,PNEUTR(5),PPHO(5,NCHARB))
      QOLD=PNEUTR(3)
      GNEUT=(QNEW**2+QOLD**2+MNESQR)/(QNEW*QOLD+SQRT((QNEW**2+MNESQR)*
     &(QOLD**2+MNESQR)))
      IF (GNEUT.LT.1.D0) THEN
        DATA=0.D0
        CALL PHOERR(4,'PHOKIN',DATA)
      ENDIF
      PARNE=GNEUT-SQRT(MAX(GNEUT**2-1.0D0,0.D0))
C--
C--   2) ...reductive boost...
      CALL PHOBO3(PARNE,PNEUTR)
C--
C--   ...calculate photon energy in the reduced system...
      NPHO=NPHO+1
      ISTPHO(NPHO)=1
      IDPHO(NPHO) =22
C--   Photon mother and daughter pointers !
      JMOPHO(1,NPHO)=IP
      JMOPHO(2,NPHO)=0
      JDAPHO(1,NPHO)=0
      JDAPHO(2,NPHO)=0
      PPHO(4,NPHO)=EPHOTO*PPHO(5,IP)/PMAVIR
C--
C--   ...and photon momenta
      CCOSTH=-COSTHG
      SSINTH=SINTHG
      TH3=PHOAN2(CCOSTH,SSINTH)
      FI3=TWOPI*PHORAN(FI3DUM)
      PPHO(1,NPHO)=PPHO(4,NPHO)*SINTHG*COS(FI3)
      PPHO(2,NPHO)=PPHO(4,NPHO)*SINTHG*SIN(FI3)
C--
C--   Minus sign because axis opposite direction of charged particle !
      PPHO(3,NPHO)=-PPHO(4,NPHO)*COSTHG
      PPHO(5,NPHO)=0.D0
C--
C--   Rotate in order to get photon along z-axis
      CALL PHORO3(-FI3,PNEUTR(1))
      CALL PHORO3(-FI3,PPHO(1,NPHO))
      CALL PHORO2(-TH3,PNEUTR(1))
      CALL PHORO2(-TH3,PPHO(1,NPHO))
      ANGLE=EPHOTO/PPHO(4,NPHO)
C--
C--   Boost to the rest frame of decaying particle
      CALL PHOBO3(ANGLE,PNEUTR(1))
      CALL PHOBO3(ANGLE,PPHO(1,NPHO))
C--
C--   Back in the parent rest frame but PNEUTR not yet oriented !
      FI4=PHOAN1(PNEUTR(1),PNEUTR(2))
      TH4=PHOAN2(PNEUTR(3),SQRT(PNEUTR(1)**2+PNEUTR(2)**2))
      CALL PHORO3(FI4,PNEUTR(1))
      CALL PHORO3(FI4,PPHO(1,NPHO))
C--
        DO 60 I=2,4
   60   PVEC(I)=0.D0
        PVEC(1)=1.D0
        CALL PHORO3(-FI3,PVEC)
        CALL PHORO2(-TH3,PVEC)
        CALL PHOBO3(ANGLE,PVEC)
        CALL PHORO3(FI4,PVEC)
        CALL PHORO2(-TH4,PNEUTR)
        CALL PHORO2(-TH4,PPHO(1,NPHO))
        CALL PHORO2(-TH4,PVEC)
        FI5=PHOAN1(PVEC(1),PVEC(2))
C--
C--   Charged particle restores original direction
        CALL PHORO3(-FI5,PNEUTR)
        CALL PHORO3(-FI5,PPHO(1,NPHO))
        CALL PHORO2(TH1,PNEUTR(1))
        CALL PHORO2(TH1,PPHO(1,NPHO))
        CALL PHORO3(FI1,PNEUTR)
        CALL PHORO3(FI1,PPHO(1,NPHO))
C--   See whether neutral system has multiplicity larger than 1...
      IF ((JDAPHO(2,IP)-JDAPHO(1,IP)).GT.1) THEN
C--   Find pointers to components of 'neutral' system
C--
        FIRST=NEUDAU
        LAST=JDAPHO(2,IP)
        DO 70 I=FIRST,LAST
          IF (I.NE.NCHARB.AND.(JMOPHO(1,I).EQ.IP)) THEN
C--
C--   Reconstruct kinematics...
            CALL PHORO3(-FI1,PPHO(1,I))
            CALL PHORO2(-TH1,PPHO(1,I))
C--
C--   ...reductive boost
            CALL PHOBO3(PARNE,PPHO(1,I))
C--
C--   Rotate in order to get photon along z-axis
            CALL PHORO3(-FI3,PPHO(1,I))
            CALL PHORO2(-TH3,PPHO(1,I))
C--
C--   Boost to the rest frame of decaying particle
            CALL PHOBO3(ANGLE,PPHO(1,I))
C--
C--   Back in the parent rest-frame but PNEUTR not yet oriented.
            CALL PHORO3(FI4,PPHO(1,I))
            CALL PHORO2(-TH4,PPHO(1,I))
C--
C--   Charged particle restores original direction
            CALL PHORO3(-FI5,PPHO(1,I))
            CALL PHORO2(TH1,PPHO(1,I))
            CALL PHORO3(FI1,PPHO(1,I))
          ENDIF
   70   CONTINUE
      ELSE
C--
C--   ...only one 'neutral' particle in addition to photon!
        DO 80 J=1,4
   80   PPHO(J,NEUDAU)=PNEUTR(J)
      ENDIF
C--
C--   All 'neutrals' treated, fill /PHOEVT/ for charged particle...
      DO 90 J=1,3
   90 PPHO(J,NCHARB)=-(PPHO(J,NPHO)+PNEUTR(J))
      PPHO(4,NCHARB)=PPHO(5,IP)-(PPHO(4,NPHO)+PNEUTR(4))
C--
      END

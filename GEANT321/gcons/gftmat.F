*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:14  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 11/01/95  08.57.12  by  S.Ravndal
*-- Author :
      SUBROUTINE GFTMAT(IMATE,IPATT,CHMECA,KDIM,TKIN,VALUE,PCUT,IXST)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       FETCH and INTERPOLATE the DE/DX and Cross sections       *
C.    *       tabulated in JMATE banks coresponding to  :              *
C.    *       material IMATE, particle IPATT, mecanism name CHMECA,    *
C.    *       kinetic energies TKIN                                    *
C.    *                                                                *
C.    *      The CHMECAnism name can be :                              *
C.    *      'HADF'   'INEF'   'ELAF'   'FISF'   'CAPF'                *
C.    *      'HADG'   'INEG'   'ELAG'   'FISG'   'CAPG'                *
C.    *      'LOSS'   'PHOT'   'ANNI'   'COMP'   'BREM'                *
C.    *      'PAIR'   'DRAY'   'PFIS'   'RAYL'   'HADG'                *
C.    *      'MUNU'   'RANG'   'STEP'                                  *
C.    *                                                                *
C.    *       For Hadronic particles it also computes the              *
C.    *       hadronic cross section from FLUKA ( '***F' ) or          *
C.    *       GHEISHA ( '***G' ) programs:                             *
C.    *       HADF or HADG -- total                                    *
C.    *       INEF or INEG -- inelastic                                *
C.    *       ELAF or ELAG -- elastic                                  *
C.    *       FISF or FISG -- fission (0.0 for FLUKA)                  *
C.    *       CAPF or CAPG -- neutron capture (0.0 for FLUKA)          *
C.    *                                                                *
C.    *             Input parameters                                   *
C.    *  IMATE   Geant material number                                 *
C.    *  IPATT   Geant particle number                                 *
C.    *  CHMECA   mechanism name of the bank to be fetched             *
C.    *  KDIM   dimension of the arrays TKIN , VALUE                   *
C.    *  TKIN   array of kinetic energy of incident particle (in Gev)  *
C.    *                                                                *
C.    *             Output parameters                                  *
C.    *  VALUE  array of energy loss (in Mev/cm) ,                     *
C.    *               or stopping range (cm) ,or continuous step (cm)  *
C.    *               or macroscopic cross section (in 1/cm)           *
C.    *  PCUT(5)  array of the physical cuts in material IMATE  (Gev)  *
C.    *  IXST   flag = 1 if the array VALUE is filled ,  =0 otherwise  *
C.    *                                                                *
C.    *    ==>Called by : <USER>  GPLMAT  GRPMAT                       *
C.    *       Authors    R.Brun, M.Maire    *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconst.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gsecti.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gctmed.inc"
#include "geant321/gfkdis.inc"
#include "geant321/gcking.inc"
#include "geant321/gckine.inc"
*
      LOGICAL CERKOV
      CHARACTER*4  CHMECA
      DIMENSION TKIN(KDIM), VALUE(KDIM), PCUT(5)
#include "geant321/gcnmec.inc"
*
*     ------------------------------------------------------------------
*
      IXST = 0
      IF(KDIM.LE.0) GO TO 999
      IMECA = 0
      DO 10  KMECA=1,NMECA
         IF(CHMECA.EQ.CHNMEC(KMECA)) THEN
            IMECA = KMECA
         ENDIF
   10 CONTINUE
      IF(IMECA.EQ.0) THEN
         WRITE(CHMAIL,'('' *** GFTMAT: Mechanism '',A,                 '
     +   //'    ''not implemented'')') CHMECA
         CALL GMAIL(0,0)
         GOTO 999
      ENDIF
      DO 20  IDIM=1, KDIM
         VALUE(IDIM)=0.
   20 CONTINUE
      DO 30  ICUT=1,5
         PCUT(ICUT)=0.
   30 CONTINUE
*
      IF(JMATE.LE.0) GO TO 999
      IF(IMATE.LE.0) GO TO 999
      IF(IMATE.GT.NMATE) GO TO 110
      JMA = LQ(JMATE-IMATE)
      IF(JMA.LE.0) GO TO 110
      A      =  Q(JMA+6)
      Z      =  Q(JMA+7)
      IF(Z.LT.1.) GO TO 999
      DENS   =  Q(JMA+8)
      RADL   =  Q(JMA+9)
      NLM    =  Q(JMA+11)
      JPROB  = LQ(JMA-4)
      AZRO   =  Q(JPROB+8)
      AHEFF  =  A
      IF(NLM .GT.1) THEN
         JMIXT = LQ(JMA-5)
         JMI1  = LQ(JMIXT-1)
         AHEFF =  Q(JMI1+1)
      ENDIF
*
      IF(JTMED.LE.0) GO TO 999
      IF(NTMED.LE.0) GO TO 999
      JBANK = JTMED
      DO 40 ITM = 1,NTMED
         JTM = LQ(JTMED-ITM)
         IF(JTM.LE.0) GO TO 40
         JTMN =  0
         IMAT =  Q(JTM+6)
         IF(IMAT.EQ.IMATE) THEN
            JTMN = LQ(JTM)
            IF(JTMN.NE.0) JBANK = JTMN
            GO TO 50
         ENDIF
   40 CONTINUE
   50 CALL UCOPY( Q(JBANK+6),PCUT(1),5)
      CUTHAD = Q(JBANK+ 4)
      ILOSS  = Q(JBANK+21)
      IMULS  = Q(JBANK+22)
      IFIELD = Q(JTM +  8)
      FIELDM = Q(JTM +  9)
      TMAXFD = Q(JTM + 10)
      STEMAX = Q(JTM + 11)
      DEEMAX = Q(JTM + 12)
      STMIN  = Q(JTM + 14)
*
      IF(JPART.LE.0) GO TO 999
      IF(IPATT.LE.0) GO TO 999
      IF(IPATT.GT.NPART) GO TO 110
      JPA = LQ(JPART-IPATT)
      IF(JPA.LE.0) GO TO 110
      ITYPE  =  Q(JPA+6)
      AMASS  =  Q(JPA+7)
      CHARGE =  Q(JPA+8)
*
* *** Find the correct pointer
*
      JBANK  = 0
      ISHIF  = 0
      RMASS  = 1.
*
* *** Photons
*
      IF (ITYPE.EQ.1) THEN
         IF (CHMECA.EQ.'PHOT') JBANK = LQ(JMA- 6)
         IF (CHMECA.EQ.'COMP') JBANK = LQ(JMA- 8)
         IF (CHMECA.EQ.'PAIR') JBANK = LQ(JMA-10)
         IF (CHMECA.EQ.'PFIS') JBANK = LQ(JMA-12)
         IF (CHMECA.EQ.'RAYL') JBANK = LQ(JMA-13)
*
* *** Electrons / positons
*
      ELSE IF (ITYPE.EQ.2) THEN
         IF (CHMECA.EQ.'LOSS') THEN
            JBANK = LQ(JMA- 1)
            IF (CHARGE.GT.0.) ISHIF = NEK1
         ELSE IF (CHMECA.EQ.'RANG') THEN
            JBANK = LQ(JMA- 15)
            IF (CHARGE.GT.0.) ISHIF = NEK1
         ELSE IF (CHMECA.EQ.'STEP') THEN
            JBANK = LQ(JTM- 1)
         ELSE IF ((CHMECA.EQ.'ANNI').AND.(CHARGE.GT.0.)) THEN
            JBANK = LQ(JMA- 7)
         ELSE IF (CHMECA.EQ.'BREM') THEN
            JBANK = LQ(JMA- 9)
         ELSE IF (CHMECA.EQ.'DRAY') THEN
            JBANK = LQ(JMA-11)
            IF (CHARGE.GT.0.) ISHIF = NEK1
         ENDIF
*
* *** Neutral hadrons (***F for FLUKA cross sections
*                      ***G for GHEISHA cross sections
*                      LOWN and N*** for MICAP cross sections)
*
      ELSE IF (ITYPE.EQ.3) THEN
         IF((CHMECA.EQ.'HADF').OR.(CHMECA.EQ.'INEF')
     +   .OR.(CHMECA.EQ.'ELAF') .OR.(CHMECA.EQ.'FISF')
     +   .OR.(CHMECA.EQ.'CAPF')) THEN
            JBANK = -3
            IF (IFINIT(5) .EQ. 0) CALL FLINIT
         ELSE IF((CHMECA.EQ.'HADG').OR.(CHMECA.EQ.'INEG')
     +   .OR.(CHMECA.EQ. 'ELAG').OR.(CHMECA.EQ.'FISG')
     +   .OR.(CHMECA.EQ.'CAPG')) THEN
            JBANK = -4
            CALL GHEINI
         ELSE IF(IMECA.GE.IBLOWN.AND.IPATT.EQ.13) THEN
            IF (IFINIT(7) .EQ. 0) CALL GMORIN
            JBANK = -5
         ENDIF
         K0OLD = K0FLAG
*
* *** Charged hadrons (***F for FLUKA cross sections
*                      ***G for GHEISHA cross sections)
* *** Heavy ions
*
      ELSE IF (ITYPE.EQ.4.OR.ITYPE.EQ.8) THEN
         RMASS = PMASS/AMASS
         IF (CHMECA.EQ.'LOSS') THEN
            JBANK = LQ(JMA- 3)
         ELSE IF (CHMECA.EQ.'RANG') THEN
            JBANK = LQ(JMA- 16) + NEK1
         ELSE IF (CHMECA.EQ.'STEP') THEN
            JBANK = -1
            JRANG = LQ(JMA -16) + NEK1
            CUTPRO = CUTHAD*RMASS
            CUTPRO = MAX(ELOW(1), MIN( CUTPRO, ELOW(NEK1)*0.99))
            IKCUT  = GEKA*LOG10(CUTPRO) + GEKB
            GKC = (CUTPRO - ELOW(IKCUT))/(ELOW(IKCUT+1) - ELOW(IKCUT))
            STOPC = (1.-GKC)*Q(JRANG+IKCUT) + GKC*Q(JRANG+IKCUT+1)
         ELSE IF (CHMECA.EQ.'DRAY') THEN
            JBANK = -2
            JPROB = LQ(JMA-4)
            AZRO  =  Q(JPROB+17)
            DCUTM =  PCUT(4)
         ELSE IF(((CHMECA.EQ.'HADF').OR.(CHMECA.EQ.'INEF')
     +   .OR.(CHMECA.EQ. 'ELAF').OR.(CHMECA.EQ.'FISF')
     +   .OR.(CHMECA.EQ.'CAPF')).AND.ITYPE.NE.8) THEN
            JBANK = -3
            IF (IFINIT(5) .EQ. 0) CALL FLINIT
         ELSE IF(((CHMECA.EQ.'HADG').OR.(CHMECA.EQ.'INEG')
     +   .OR.(CHMECA.EQ. 'ELAG').OR.(CHMECA.EQ.'FISG')
     +   .OR.(CHMECA.EQ.'CAPG')).AND.ITYPE.NE.8) THEN
            JBANK = -4
            CALL GHEINI
         ENDIF
         K0OLD = K0FLAG
*
* *** Muons
*
      ELSE IF (ITYPE.EQ.5) THEN
         IF (CHMECA.EQ.'LOSS') THEN
            JBANK = LQ(JMA- 2)
         ELSE IF (CHMECA.EQ.'RANG') THEN
            JBANK = LQ(JMA- 16)
         ELSE IF (CHMECA.EQ.'STEP') THEN
            JBANK = LQ(JTM- 2)
         ELSE IF (CHMECA.EQ.'MUNU') THEN
            JBANK = LQ(JMA- 14)
         ELSE IF (CHMECA.EQ.'BREM') THEN
            JBANK = LQ(JMA- 9)
            ISHIF = 2*NEK1
         ELSE IF (CHMECA.EQ.'PAIR') THEN
            JBANK = LQ(JMA- 10)
            ISHIF = NEK1
         ELSE IF (CHMECA.EQ.'DRAY') THEN
            JBANK = LQ(JMA-11)
            ISHIF = 2*NEK1
         ENDIF
*
* *** Geantinos
*
      ELSEIF (ITYPE.EQ.6) THEN
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
         JBANK = 0
*
* *** Cerenkov
*
      ELSEIF (ITYPE.EQ.7) THEN
         IF (CHMECA.EQ.'LABS') THEN
*
* *** Not implemented yet!
            JBANK=0
         ENDIF
*
      ENDIF
      CERKOV=.FALSE.
      IF(CHARGE.NE.0.AND.ITCKOV.NE.0) THEN
         IF(IQ(JTM-2).GE.3) THEN
            IF(LQ(JTM-3).NE.0.AND.LQ(LQ(JTM-3)-3).NE.0) THEN
*
* ***  In this tracking medium Cerenkov photons are generated and
* ***  tracked. Set to 1 the corresponding flag and calculate the
* ***  relevant pointers.
*
               CERKOV = .TRUE.
               JTCKOV = LQ(JTM-3)
               JABSCO = LQ(JTCKOV-1)
               JEFFIC = LQ(JTCKOV-2)
               JINDEX = LQ(JTCKOV-3)
               JCURIN = LQ(JTCKOV-4)
               NPCKOV = Q(JTCKOV+1)
            ENDIF
         ENDIF
      ENDIF
 
      IF(JBANK.EQ.0) GO TO 999
      IXST = 1
*
*
      JBANK = JBANK + ISHIF
      DO 100 IKB = 1,KDIM
*        Find bin number in table JMATE
         EKP = TKIN(IKB)*RMASS
         EKP = MAX(ELOW(1), MIN( EKP, ELOW(NEK1)*0.99))
         IKP=GEKA*LOG10(EKP)+GEKB +0.001
         GKRA=(EKP-ELOW(IKP))/(ELOW(IKP+1)-ELOW(IKP))
*
         IF(JBANK.GT.0) THEN
*           Retieve value from bank JMATE
            VALUE(IKB) =  (1.-GKRA)*Q(JBANK+IKP) + GKRA*Q(JBANK+IKP+1)
            IF ((CHMECA.EQ.'PHOT').AND.(EKP.GE.0.05 )) VALUE(IKB) =
     +      BIG
            IF ((CHMECA.EQ.'MUNU').AND.(EKP.LT.0.05 )) VALUE(IKB) =
     +      BIG
            IF ( CHMECA.EQ.'LOSS') THEN
               VALUE(IKB) = VALUE(IKB)*CHARGE**2*1.E+3
            ELSE IF (CHMECA.EQ.'RANG') THEN
               VALUE(IKB) = VALUE(IKB)/(RMASS*CHARGE*CHARGE)
            ELSE IF (CHMECA.NE.'STEP') THEN
               IF (VALUE(IKB).GT.0.) THEN
                  VALUE(IKB) = 1./VALUE(IKB)
               ELSE
                  VALUE(IKB) = 1./BIG
               ENDIF
            ENDIF
*
         ELSEIF (JBANK.EQ.-1) THEN
*           Compute step due to muls + loss + field
            GEKIN=TKIN(IKB)
            GETOT=GEKIN+AMASS
            PMOM =SQRT(GEKIN*(GETOT+AMASS))
            SFIELD = BIG
            SMULS  = BIG
            SLOSS  = BIG
            STOPMX = BIG
            IF (IFIELD*FIELDM.NE.0.)
     +         SFIELD = 3333.*DEGRAD*TMAXFD*PMOM/ABS(FIELDM*CHARGE)
            IF (IMULS.GT.0.)
     +         SMULS = MIN (2232.*RADL*((PMOM**2)/(GETOT*CHARGE))**2 ,
     +                      10.*RADL )
            IF (ILOSS*DEEMAX.GT.0.) THEN
               STOPP  = (1.-GKRA)*Q(JRANG+IKP) + GKRA*Q(JRANG+IKP+1)
               STOPMX = (STOPP - STOPC)/(RMASS*CHARGE*CHARGE)
               IF (STOPMX.LT.0.) STOPMX = 0.
               EKF = MAX ( ELOW(1) , (1.-DEEMAX)*EKP )
               IKF = GEKA*LOG10(EKF) + GEKB
               GKF = (EKF-ELOW(IKF))/(ELOW(IKF+1)-ELOW(IKF))
               SLOSP= STOPP - (1.-GKF)*Q(JRANG+IKF) - GKF*Q(JRANG+IKF+1)
               SLOSS = SLOSP/(RMASS*CHARGE*CHARGE)
            ENDIF
            IF(CERKOV) THEN
               VECT(7)=SQRT(TKIN(IKB)*(TKIN(IKB)+2*AMASS))
               CALL GNCKOV
               STCKOV = MXPHOT/MAX(3.*DNDL,1E-10)
            ELSE
               STCKOV=BIG
            ENDIF
*
            IF (STOPMX.LE.STMIN) THEN
               VALUE(IKB) = STOPMX
            ELSE
               VALUE(IKB) = MAX(STMIN,MIN(STCKOV,SLOSS,SFIELD,SMULS,
     +         STEMAX))
            ENDIF
*
         ELSEIF (JBANK.EQ.-2) THEN
*           Compute delta ray cross section for hadrons
            GEKIN=TKIN(IKB)
            GETOT=GEKIN+AMASS
            GAMASS=GETOT+AMASS
            BET2=GEKIN*GAMASS/(GETOT*GETOT)
            TMAX=EMASS*GEKIN*GAMASS/(0.5*AMASS*AMASS+EMASS*GETOT)
            IF(TMAX.GT.DCUTM)THEN
               Y=DCUTM/TMAX
               SIG=(1.-Y+BET2*Y*LOG(Y))/DCUTM
               IF(AMASS.GT.0.9)SIG=SIG+0.5*(TMAX-DCUTM)/(GETOT*GETOT)
               VALUE(IKB)=SIG*AZRO*CHARGE*CHARGE*EMASS/BET2
            ELSE
               VALUE(IKB)=1./BIG
            ENDIF
*
*           compute hadronic cross section from FLUKA code
*
         ELSEIF (JBANK.EQ.-3) THEN
            GEKIN=TKIN(IKB)
            PMOM = SQRT(GEKIN*(GEKIN+2*AMASS))
            NMAT = IMATE
            VECT(7) = PMOM
            IF (IPATT.NE.IPART) THEN
                IOLDP = IPART
                IPART = IPATT
                CALL FLDIST
                IPART = IOLDP
            ELSE
                CALL FLDIST
            END IF
            IF (CHMECA.EQ.'HADF') VALUE(IKB)= FSIG
            IF (CHMECA.EQ.'INEF') VALUE(IKB)= SINE
            IF (CHMECA.EQ.'ELAF') VALUE(IKB)= SELA
            IF (CHMECA.EQ.'FISF') VALUE(IKB)= 0.0
            IF (CHMECA.EQ.'CAPF') VALUE(IKB)= 0.0
*
*           compute hadronic cross section from GHEISHA code
*
         ELSEIF (JBANK.EQ.-4) THEN
            GEKIN=TKIN(IKB)
            PMOM = SQRT(GEKIN*(GEKIN+2*AMASS))
            K0FLAG = 1
*           (compounds)
            IF(NLM.GT.1) THEN
               HHHH=0.
               IF (JTMN.GT.0) HHHH=Q(JTMN+26)
               VALUE(IKB)=GHESIG(PMOM,GEKIN,A,Q(JMIXT+1), Q(JMIXT+NLM+
     +         1),Q(JMIXT+2*NLM+1),NLM,DENS,HHHH,IPATT)
               IF (CHMECA .EQ. 'INEG') THEN
                  VALUE(IKB) = 0.
                  DO 60 K=1,NLM
                     VALUE(IKB) = AIIN(K) + VALUE(IKB)
   60             CONTINUE
               ELSE IF (CHMECA .EQ. 'ELAG') THEN
                  VALUE(IKB) = 0.
                  DO 70 K=1,NLM
                     VALUE(IKB) = AIEL(K) + VALUE(IKB)
   70             CONTINUE
               ELSE IF (CHMECA .EQ. 'FISG') THEN
                  VALUE(IKB) = 0.
                  DO 80 K=1,NLM
                     VALUE(IKB) = AIFI(K) + VALUE(IKB)
   80             CONTINUE
               ELSE IF (CHMECA .EQ. 'CAPG') THEN
                  VALUE(IKB) = 0.
                  DO 90 K=1,NLM
                     VALUE(IKB) = AICA(K) + VALUE(IKB)
   90             CONTINUE
               ENDIF
            ELSE
*           (simple elements)
               VALUE(IKB)=GHESIG(PMOM,GEKIN,A,A,Z,1.,1,DENS,0.,IPATT)
               IF (CHMECA .EQ. 'INEG') VALUE(IKB) = AIIN(1)
               IF (CHMECA .EQ. 'ELAG') VALUE(IKB) = AIEL(1)
               IF (CHMECA .EQ. 'FISG') VALUE(IKB) = AIFI(1)
               IF (CHMECA .EQ. 'CAPG') VALUE(IKB) = AICA(1)
            END IF
            K0FLAG = K0OLD
*
*           compute the cross-section for low-energy neutrons
*           from MICAP code
*
         ELSEIF (JBANK.EQ.-5) THEN
            GEKIN=TKIN(IKB)
            IF (GEKIN.LE..02) THEN
               IF (CHMECA .EQ. 'LOWN') THEN
                  VALUE(IKB) = SIGMOR(GEKIN*1.E+9,IMATE)
               ELSE
                  NMAT = IMATE
                  CALL GMXSEC (IMECA,VALUE(IKB))
               ENDIF
            ELSE
               VALUE(IKB) = 0.
            ENDIF
*
         ENDIF
  100 CONTINUE
*
      GO TO 999
  110 WRITE(CHMAIL,10100) IMATE ,IPATT
      CALL GMAIL(0,0)
*
10000 FORMAT(' ***** GFTMAT: No processes active for geantinos')
10100 FORMAT(' ***** GFTMAT error : material',I4,
     +       '  or particle',I4,' not defined'   )
  999 END

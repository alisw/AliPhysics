*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 13/12/94  17.08.38  by  S.Giani
*-- Author :
      SUBROUTINE GPRMAT(IMATE,IPART,MECAN,KDIN,TKIN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       INTERPOLATE and PRINT the DE/DX ,stopping range and      *
C.    *       Cross sections tabulated in JMATE banks corresponding to *
C.    *       material IMATE, particle IPART, mecanism name MECAN ,    *
C.    *       kinetic energies TKIN.                                   *
C.    *                                                                *
C.    *      The MECAnism name can be :                                *
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
C.    *  IPART   Geant particle number                                 *
C.    *  MECAN   mechanism name of the bank to be fetched              *
C.    *  KDIM   dimension of the arrays TKIN , VALUE                   *
C.    *  TKIN   array of kinetic energy of incident particle (in Gev)  *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Authors    R.Brun, M.Maire    *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      PARAMETER (MMX=100)
      CHARACTER*(*) MECAN
      CHARACTER*4  MECA
      CHARACTER*4  KU1 , KU2 , KU3 , KU(5)
      DIMENSION   TKIN(KDIN),VALUE(MMX),SIGT(MMX),PCUT(5),CU(5)
*
#include "geant321/gcnmec.inc"
*
*     ------------------------------------------------------------------
*
      KDIM = MIN(KDIN,MMX)
      IF (KDIM.LE.0) GO TO 999
*
      IF (JMATE.LE.0) GO TO 999
      IF (IMATE.LE.0) GO TO 999
      IF (IMATE.GT.NMATE) GO TO 90
      JMA = LQ(JMATE-IMATE)
      IF (JMA.LE.0) GO TO 90
*
      IF (JPART.LE.0) GO TO 999
      IF (IPART.LE.0) GO TO 999
      IF (IPART.GT.NPART) GO TO 90
      JPA = LQ(JPART-IPART)
      IF (JPA.LE.0) GO TO 90
*
      DO 10 JSIG=1,MMX
         SIGT(JSIG)=0.
   10 CONTINUE
      IF(MECAN.EQ.'ALL') THEN
         N1 = 1
         N2 = NMECA
      ELSE
         N1 = 0
         DO 20  IMECA=1,NMECA
            IF(MECAN.EQ.CHNMEC(IMECA)) THEN
               N1 = IMECA
            ENDIF
   20    CONTINUE
         IF(N1.EQ.0) THEN
            WRITE(CHMAIL,'('' *** GPRMAT: Mechanism '',A,
     +      '' not implemented'')') MECAN
            CALL GMAIL(0,0)
            GOTO 999
         ENDIF
         N2 = N1
      ENDIF
      DO 60  IMEC = N1,N2
         IF(CHNMEC(IMEC).NE.'NULL') THEN
            MECA = CHNMEC(IMEC)
            CALL GFTMAT(IMATE,IPART,MECA,KDIM,TKIN,VALUE,PCUT,IXST)
            IF(IXST.EQ.0) GO TO 60
            CHMAIL='1'
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10100)(IQ(JMA+I),I=1,5),MECA,(IQ(JPA+J),J=1,5)
            CALL GMAIL(0,0)
            CHMAIL(31:)='-----------------------------------------'
            CALL GMAIL(0,1)
            CHMAIL=' '
            DO 30  K=1,5
   30       CALL GEVKEV(PCUT(K),CU(K),KU(K))
            WRITE(CHMAIL,10200) (CU(K),KU(K),K=1,5)
            CALL GMAIL(0,1)
*
            IF (MECA.EQ.'LOSS'.OR.MECA.EQ.'RANG'.OR.MECA.EQ.'STEP')
     +      THEN
               IF (MECA.EQ.'LOSS') WRITE(CHMAIL,10300)
               IF (MECA.EQ.'RANG') WRITE(CHMAIL,10400)
               IF (MECA.EQ.'STEP') WRITE(CHMAIL,10500)
               CALL GMAIL(0,1)
               NROW = (KDIM-1)/3 + 1
               DO 40  IKB=1,NROW
                  IK = IKB
                  DE1 = VALUE(IK)
                  CALL GEVKEV(TKIN(IK),EK1,KU1)
*
                  IK = IKB + NROW
                  IF (IK.GT.KDIM) IK=KDIM
                  DE2 = VALUE(IK)
                  CALL GEVKEV(TKIN(IK),EK2,KU2)
*
                  IK = IKB + 2*NROW
                  IF (IK.GT.KDIM) IK=KDIM
                  DE3 = VALUE(IK)
                  CALL GEVKEV(TKIN(IK),EK3,KU3)
*
                  WRITE(CHMAIL,10600) EK1,KU1,DE1,EK2,KU2,DE2,EK3,KU3,
     +            DE3
                  CALL GMAIL(0,0)
   40          CONTINUE
            ELSE
               WRITE(CHMAIL,10700)
               CALL GMAIL(0,1)
               NROW = (KDIM-1)/2 + 1
               DO 50  IKB=1,NROW
                  IK = IKB
                  SIG1 = VALUE(IK)
                  AL1=0.
                  IF(SIG1.NE.0.)AL1 = 1./SIG1
                  SIGT(IK) = SIGT(IK) + SIG1
                  CALL GEVKEV(TKIN(IK),EK1,KU1)
*
                  IK = IKB + NROW
                  IF (IK.GT.KDIM) IK=KDIM
                  SIG2 = VALUE(IK)
                  AL2=0.
                  IF(SIG2.NE.0.)AL2 = 1./SIG2
                  SIGT(IK) = SIGT(IK) + SIG2
                  CALL GEVKEV(TKIN(IK),EK2,KU2)
*
                  WRITE(CHMAIL,10900) EK1,KU1,SIG1,AL1,EK2,KU2,SIG2,AL2
                  CALL GMAIL(0,0)
   50          CONTINUE
            ENDIF
         ENDIF
   60 CONTINUE
*
* *** print total cross section
      IF (MECAN.EQ.'ALL') THEN
         MECA = 'SIGT'
         CHMAIL='1'
         CALL GMAIL(0,0)
         WRITE(CHMAIL,10100)(IQ(JMA+I),I=1,5),MECA,(IQ(JPA+J),J=1,5)
         CALL GMAIL(0,0)
         CHMAIL(31:)='-----------------------------------------'
         CALL GMAIL(0,1)
         CHMAIL=' '
         DO 70  K=1,5
   70    CALL GEVKEV(PCUT(K),CU(K),KU(K))
         WRITE(CHMAIL,10200) (CU(K),KU(K),K=1,5)
         CALL GMAIL(0,1)
         WRITE(CHMAIL,10800)
         CALL GMAIL(0,1)
         NROW = (KDIM-1)/2 + 1
         DO 80  IKB=1,NROW
            IK = IKB
            SIG1 = SIGT(IK)
            AL1=0.
            IF(SIG1.NE.0.)AL1 = 1./SIG1
            CALL GEVKEV(TKIN(IK),EK1,KU1)
*
            IK = IKB + NROW
            IF (IK.GT.KDIM) IK=KDIM
            SIG2 = SIGT(IK)
            AL2=0.
            IF(SIG2.NE.0.)AL2 = 1./SIG2
            CALL GEVKEV(TKIN(IK),EK2,KU2)
*
            WRITE(CHMAIL,10900) EK1,KU1,SIG1,AL1,EK2,KU2,SIG2,AL2
            CALL GMAIL(0,0)
   80    CONTINUE
      ENDIF
*
      GO TO 999
*
   90 WRITE(CHMAIL,10000) IMATE ,IPART
      CALL GMAIL(0,0)
*
10000 FORMAT(' ***** GPRMAT error : material',I4,
     +       '  or particle',I4,' not defined'   )
10100 FORMAT(30X,5A4,A4, ' for  ',5A4)
10200 FORMAT(  6X,'BCUTE =',F6.2,A4,3X,'BCUTM =',F6.2,A4,3X,
     +             'DCUTE =',F6.2,A4,3X,'DCUTM =',F6.2,A4,3X,
     +            'PPCUTM =',F6.2,A4 )
10300 FORMAT(  6X,'kinetic energy   DE/DX(mev/cm)',
     +         6X,'kinetic energy   DE/DX(mev/cm)',
     +         6X,'kinetic energy   DE/DX(mev/cm)')
10400 FORMAT(  6X,'kinetic energy   Stop range cm',
     +         6X,'kinetic energy   Stop ramge cm',
     +         6X,'kinetic energy   Stop range cm')
10500 FORMAT(  6X,'kinetic energy   Mulof step cm',
     +         6X,'kinetic energy   Mulof step cm',
     +         6X,'kinetic energy   Mulof step cm')
10600 FORMAT( 3(F16.2,A4,E15.4))
10700 FORMAT(  6X,'kinetic energy   Sigma (1/cm)    Lambda (cm)',
     +         6X,'kinetic energy   Sigma (1/cm)    Lambda (cm)')
10800 FORMAT(  6X,'kinetic energy   Sigto (1/cm)    Lambda (cm)',
     +         6X,'kinetic energy   Sigto (1/cm)    Lambda (cm)')
10900 FORMAT( 2(F16.2,A4,2(E15.4)))
  999 END

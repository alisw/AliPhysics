*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 10/01/95  17.36.36  by  S.Ravndal
*-- Author :
      SUBROUTINE GPLMAT(IMATES,IPART,MECAN,KDIN,TKIN,IDM)
C.
C     ******************************************************************
C.    *                                                                *
C.    *       INTERPOLATE and PLOT  the DE/DX and Cross sections       *
C.    *       tabulated in JMATE banks corresponding to :              *
C.    *       material IMATE, particle IPART, mecanism name HMECAN,    *
C.    *       kinetic energies TKIN                                    *
C.    *                                                                *
C.    *      The MECANism name can be :                                *
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
C.    *  MECAN  mechanism name of the bank to be fetched               *
C.    *  KDIM   dimension of the arrays TKIN , VALUE                   *
C.    *  TKIN   array of kinetic energy of incident particle (in Gev)  *
C.    *  IDM    convention for histogramming mode :                    *
C.    *         IDM.gt.0  fill , print , keep histogram(s)             *
C.    *         IDM.eq.0  fill , print , delete histogram(s)           *
C.    *         IDM.lt.0  fill , noprint , keep histogram(s)           *
C.    *           The histogram IDentificator will be :                *
C.    *             10000*imate + 100*ipart + imeca                    *
C.    *          where IMECA is the link number in stucture JMATE      *
C.    *          (see Geant3 writeup CONS 199)                         *
C.    *           for 'HADG'  imeca = 17                               *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Authors    R.Brun, M.Maire    *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcphys.inc"
      PARAMETER (MMX= 201,NCOL= 5)
      CHARACTER*(*) MECAN
      CHARACTER*4 MECA , KU(NCOL)
      CHARACTER   NAPART*16 , NAMATE*16 , CHTITL*68
      DIMENSION   TKIN(KDIN),VALUE(MMX),SIGT(MMX),PCUT(5)
      DIMENSION   KI(NCOL),EK(NCOL)
      LOGICAL     LXBARN
*
#include "geant321/gcnmec.inc"
*
*     ------------------------------------------------------------------
*
      IF (KDIN.LE.0)  GO TO 999
      KDIM = MIN(KDIN,MMX)
      IF(IMATES.LT.0) THEN
         LXBARN=.TRUE.
      ELSEIF(IMATES.GT.0) THEN
         LXBARN=.FALSE.
      ELSE
         GOTO 999
      ENDIF
      IMATE=ABS(IMATES)
*
      IF (JMATE.LE.0) GO TO 999
      IF (IMATE.GT.NMATE) GO TO 80
      JMA = LQ(JMATE-IMATE)
      IF  (JMA.LE.0) GO TO 80
      CALL UHTOC(IQ(JMA+1),4,NAMATE,16)
      IF(LXBARN) THEN
         CMIBAR=Q(JMA+6)/(AVO*Q(JMA+8))
      ELSE
         CMIBAR=1.
      ENDIF
*
      IF (JPART.LE.0) GO TO 999
      IF (IPART.LE.0) GO TO 999
      IF (IPART.GT.NPART) GO TO 80
      JPA = LQ(JPART-IPART)
      IF (JPA.LE.0) GO TO 80
      CALL UHTOC(IQ(JPA+1),4,NAPART,16)
*
* *** Print  bin meaning
      IF (IDM.GE.0) THEN
         CHMAIL='1'
         CALL GMAIL(0,0)
         CHMAIL=' '
         CHMAIL(31:)='Kinetic energy bin meaning'
         CALL GMAIL(0,0)
         CHMAIL(31:)='--------------------------'
         CALL GMAIL(0,1)
         NROW = (KDIM-1)/NCOL + 1
         DO 20  IR=1,NROW
            DO 10  IC=1,NCOL
               IKB = IR + (IC-1)*NROW
               IF (IKB.GT.KDIM) IKB=KDIM
               KI(IC) = IKB
               CALL GEVKEV(TKIN(IKB),EK(IC),KU(IC))
   10       CONTINUE
            WRITE(CHMAIL,10200) (KI(IC),EK(IC),KU(IC),IC=1,NCOL)
            CALL GMAIL(0,0)
   20    CONTINUE
      ENDIF
*
      BIGINV= 1000./BIG
      DO 30  JMX = 1, MMX
         SIGT(JMX) = 0.
   30 CONTINUE
      IF(MECAN.EQ.'ALL') THEN
         N1 = 1
         N2 = NMECA
      ELSE
         N1 = 0
         DO 40  IMECA=1,NMECA
            IF(MECAN.EQ.CHNMEC(IMECA)) THEN
               N1 = IMECA
            ENDIF
   40    CONTINUE
         IF(N1.EQ.0) THEN
            WRITE(CHMAIL,'('' *** GPLMAT: Mechanism '',A,
     +      '' not implemented'')') MECAN
            CALL GMAIL(0,0)
            GOTO 999
         ENDIF
         N2 = N1
      ENDIF
      DO 60  IMEC = N1,N2
C
         IF (MECAN.EQ.'ALL') THEN
             IF (CHNMEC(IMEC).EQ.'RANG') GO TO 60
             IF (CHNMEC(IMEC).EQ.'STEP') GO TO 60
         END IF
C
         IF(CHNMEC(IMEC).NE.'NULL') THEN
            MECA = CHNMEC(IMEC)
            CALL GFTMAT(IMATE,IPART,MECA,KDIM,TKIN,VALUE,PCUT,IXST)
            IF(IXST.EQ.0) GO TO 60
*
* ***    Book histogram
            ISIG = 0
            IF (MECA.EQ.'LOSS') THEN
               CHTITL = NAPART//' in '//NAMATE//'   dE/dx (MeV/cm)'
            ELSEIF (MECA.EQ.'RANG') THEN
               CHTITL = NAPART//' in '//NAMATE//'   Stopping range (cm)'
            ELSEIF (MECA.EQ.'STEP') THEN
               CHTITL = NAPART//' in '//NAMATE//'   continuous step '
     +         //'(cm)'
            ELSE
               CHTITL = NAPART//' in '//NAMATE//'   '//MECA// ' cross '
     +         //'section'
               IF(LXBARN) THEN
                  CHTITL(LNBLNK(CHTITL)+1:) = ' (barn)'
               ELSE
                  CHTITL(LNBLNK(CHTITL)+1:) = ' (1/cm)'
               ENDIF
               ISIG = 1
            ENDIF
*
            ID = 10000*IMATE + 100*IPART + IMEC
            CALL HBOOKB(ID,CHTITL,KDIM-1,TKIN,0.)
*
* ***    Fill histogram
*
            VALMI = MAX (BIGINV,VMAX(VALUE,KDIM)*1.E-8)
            DO 50  IKB = 1,KDIM
               IF (MECA.NE.'LOSS'.AND.MECA.NE.'RANG'
     +             .AND.MECA.NE.'STEP')
     +             VALUE(IKB)=VALUE(IKB)*CMIBAR
               IF (VALUE(IKB).GE.VALMI) THEN
                  CALL HFILL(ID,TKIN(IKB),0.,VALUE(IKB))
               ENDIF
               IF (ISIG.EQ.1) THEN
                  IF(MECA(1:3).NE.'INE'.AND.MECA(1:3).NE.'ELA'.AND.
     +            MECA(1:3).NE.'FIS'.AND.MECA(1:3).NE.'CAP'.AND.
     +            MECA(1:3).NE.'HAD'.AND.IMEC.LT.IBLOWN) THEN
                     SIGT(IKB) = SIGT(IKB) + VALUE(IKB)
                  ELSE IF (MECA(1:3).EQ.'HAD') THEN
                     IF ((MECA.EQ.'HADG'.AND.IHADR.LE.2).OR. (MECA.EQ.
     +               'HADF'.AND.IHADR.EQ.4)) THEN
                        SIGT(IKB) = SIGT(IKB) + VALUE(IKB)
                     ENDIF
                  ENDIF
               ENDIF
   50       CONTINUE
            CALL HIDOPT(ID,'LOGY')
            IF(IDM.GE.0) CALL HPHIST(ID,' ',0)
            IF(IDM.EQ.0) CALL HDELET(ID)
         ENDIF
   60 CONTINUE
*
* *** plot total cross section and mean free path
      IF (MECAN.EQ.'ALL') THEN
         CHTITL= NAPART//' in '//NAMATE//'   total cross section'
         IF(LXBARN) THEN
            CHTITL(LNBLNK(CHTITL)+1:) = ' (barn)'
         ELSE
            CHTITL(LNBLNK(CHTITL)+1:) = ' (1/cm)'
         ENDIF
         ID = 10000*IMATE + 100*IPART + NMECA+1
         CALL HBOOKB(ID,CHTITL,KDIM-1,TKIN,0.)
*
         CHTITL= NAPART//' in '//NAMATE//'   total mean free path (cm)'
         II = ID + 1
         CALL HBOOKB(II,CHTITL,KDIM-1,TKIN,0.)
*
         VALMI = MAX (BIGINV,VMAX( SIGT,KDIM)*1.E-8)
         DO 70  IKB = 1,KDIM
            IF (SIGT(IKB).GE.VALMI) THEN
               CALL HFILL(ID,TKIN(IKB),0.,       SIGT(IKB))
               CALL HFILL(II,TKIN(IKB),0.,CMIBAR/SIGT(IKB))
            ENDIF
   70    CONTINUE
         CALL HIDOPT(ID,'LOGY')
         IF(IDM.GE.0) CALL HPHIST(ID,' ',0)
         IF(IDM.EQ.0) CALL HDELET(ID)
*
         CALL HIDOPT(II,'LOGY')
         IF(IDM.GE.0) CALL HPHIST(II,' ',0)
         IF(IDM.EQ.0) CALL HDELET(II)
      ENDIF
*
      GO TO 999
*
   80 WRITE(CHMAIL,10000) IMATE ,IPART
      CALL GMAIL(0,0)
10000 FORMAT(' ***** GPLMAT error : material',I4,
     +       '  or particle',I4,' not defined'   )
10100 FORMAT(6X,'BCUTE =',F6.2,A4,3X,'BCUTM =',F6.2,A4,3X,
     +             'DCUTE =',F6.2,A4,3X,'DCUTM =',F6.2,A4,3X,
     +            'PPCUTM =',F6.2,A4 )
10200 FORMAT(1X,5('   bin ',I3,' =',F7.2,A4))
  999 END

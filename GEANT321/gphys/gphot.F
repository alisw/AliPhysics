*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GPHOT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  GENERATES PHOTO ELECTRIC MECHANISM                            *
C.    *  Corrected version of L. Urban's routine.                      *
C.    *  Improvements:                                                 *
C.    *    1. Angular distributions of photoelectrons from K-L3 shells *
C.    *    2. Generation of shell decay mode                           *
C.    *    3. Probability of interactioon with a shell = function      *
C.    *       of photon energy                                         *
C.    *                                                                *
C.    *    ==>CALLED BY : GTGAMA                                       *
C.    *       AUTHOR    : J. Chwastowski                               *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcking.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcunit.inc"
      DIMENSION POT(4),PROB(4),RNA(9)
      EQUIVALENCE (RNA(1),RN01),(RNA(2),RN02),(RNA(3),RN03)
      EQUIVALENCE (RNA(4),RN04),(RNA(5),RN05),(RNA(6),RN06)
      EQUIVALENCE (RNA(7),RN07),(RNA(8),RN08),(RNA(9),RN09)
      EQUIVALENCE (POT(1),POTK),(POT(2),POTL1)
      EQUIVALENCE (POT(3),POTL2),(POT(4),POTL3)
      EQUIVALENCE (PROB(1),PROBK),(PROB(2),PROBL1)
      EQUIVALENCE (PROB(3),PROBL2),(PROB(4),PROBL3)
      SAVE ZINOLD,POT,NSHELL
      DATA ZINOLD / 0.0 /
C.
C.    ------------------------------------------------------------------
C.
      KCASE = NAMEC(8)
C
C             STOP ELECTRON ?
C
C Check if the photoelectric effect was activated. If not deposit
C gamma & return
      IF(IPHOT.NE.1) THEN
         ISTOP = 2
         NGKINE= 0
         DESTEP = DESTEP + VECT(7)
         VECT(7) = 0.
         GEKIN = 0.
         GETOT = 0.
      ELSE
         E=VECT(7)
         CALL GRNDM(RNA,9)
         JPHXS = LQ(JPHOT-1)
         NZ = Q(JPHXS+1)
         IF(NZ.GT.1) THEN
            QS = 0.0
            QS2 = GPHSG1(E)*RN01
            DO 10 I = 1,NZ-1
               QS1 = GPHSGP(I,E)
               QS = QS+QS1
               IF(QS2.LE.QS) THEN
                  K = I
                  GO TO 20
               ENDIF
   10       CONTINUE
            K = NZ
   20       CONTINUE
            JPHFN = LQ(JPHXS-K)
            NUSED = Q(JPHFN+1)*5+1
            JFN = JPHFN+NUSED
            ZINT = Q(JPHXS+1+K)
         ELSE
            JPHFN = LQ(JPHXS-1)
            NUSED = Q(JPHFN+1)*5+1
            JFN = JPHFN+NUSED
            ZINT = Q(JPHXS+1+1)
         ENDIF
C COPY SHELLS POTENTIALS FROM THE ZEBRA STUCTURE
C Check if this atom was used in last entry
         IF(ZINT.NE.ZINOLD) THEN
            NSHELL = Q(JFN+1)
            DO 30 I = 1,NSHELL
               POT(I) = Q(JFN+1+I)
   30       CONTINUE
            ZINOLD = ZINT
         ENDIF
C Check if E-gamma is bigger than the L3 ionization potential.
C This will make GPHOT a little faster.
         ISHELL = 0
         PROB(1) = 0.
         PROB(2) = 0.
         PROB(3) = 0.
         PROB(4) = 0.
         IF(E.GE.POTL3) THEN
C If ZINT < 5 we can have K shell only, so
            IF(ZINT.LT.5) THEN
               IF(E.GT.POTK) THEN
                  PROBK = 1.
                  TK = E-POTK
                  ISHELL = 1
               ENDIF
            ELSE
C The probabilities given below come from crude approximation
C It uses the jump ratios and assumes that they are valid for the whole energy
C range.
               IF(E.LT.POTL2) THEN
                  PROBL3 = 1.0
                  TK = E-POTL3
                  ISHELL = 4
               ELSE
                  E3 = E-POTL3
                  GAMAL3 = E3/EMASS+1.
                  BETAL3 = SQRT(E3*(E3+2.0*EMASS))/(E+EMASS)
                  E2 = E-POTL2
                  GAMAL2 = E2/EMASS+1.
                  BETAL2 = SQRT(E2*(E2+2.0*EMASS))/(E+EMASS)
                  EFRAC = EMASS/E
                  PROBL3 = GAVRL3(GAMAL3,BETAL3,EFRAC)
                  PROBL2 = GAVRL2(GAMAL2,BETAL2,EFRAC)
                  ANOR = 1./(PROBL3+PROBL2)
                  PROBL3 = PROBL3*ANOR
                  PROBL2 = PROBL2*ANOR
                  IF(E.LT.POTL1) THEN
                     IF(RN02.LT.PROBL3) THEN
                        ISHELL = 4
                        TK = E-POTL3
                     ELSE
                        ISHELL = 3
                        TK = E-POTL2
                     ENDIF
                  ELSE
C Parametrization of L1 jump ratio gives constant 1.2
                     PROBL1 = 1.-1./1.2
                     IF(E.LT.POTK) THEN
                        PROBL2 = (1.-PROBL1)*PROBL2
                        PROBL3 = (1.-PROBL1)*PROBL3
                     ELSE
                        PROBK = 125./ZINT+3.5
                        PROBK = 1.-1/PROBK
                        PROBL1 = (1.-PROBK)*PROBL1
                        PROBL2 = (1.-PROBK-PROBL1)*PROBL2
                        PROBL3 = (1.-PROBK-PROBL1)*PROBL3
                     ENDIF
                     IF(POTL3.LE.0.0) PROBL3 = 0.0
                     IF(POTL2.LE.0.0) PROBL2 = 0.0
                     IF(POTL1.LE.0.0) PROBL1 = 0.0
                     ANOR = PROBK+PROBL1+PROBL2+PROBL3
                     IF(ANOR.GT.0.0) THEN
                        ANOR = 1./ANOR
                        PROBK = PROBK*ANOR
                        PROBL1 = PROBL1*ANOR+PROBK
                        PROBL2 = PROBL2*ANOR+PROBL1
                        PROBL3 = PROBL3*ANOR+PROBL2
                        ISHELL = 4
                        TK = E-POTL3
                        IF(RN02.LE.PROBK) THEN
                           ISHELL = 1
                           TK = E-POTK
                        ELSEIF(RN02.LE.PROBL1) THEN
                           ISHELL = 2
                           TK = E-POTL1
                        ELSEIF(RN02.LE.PROBL2) THEN
                           ISHELL = 3
                           TK = E-POTL2
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF(TK.LE.CUTELE) ISHELL = -ISHELL
         ENDIF
         IF(ISHELL.LT.1) THEN
C None of the shells was chosen because of the CUTELE
            ISTOP = 2
            IF(ISHELL.LT.0) THEN
               DESTEP = DESTEP+TK
            ELSEIF(ISHELL.EQ.0) THEN
               DESTEP = DESTEP+VECT(7)
            ENDIF
            NGKINE= 0
            VECT(7) = 0.
            GEKIN = 0.
            GETOT = 0.
         ELSE
C
C             ENERGY AND MOMENTUM OF PHOTOELECTRON
C
            EEL=TK + EMASS
            PEL=SQRT((TK+2.*EMASS)*TK)
            BETA = PEL/EEL
            ISTOP = 1
            NGKINE = 1
            IF(ISHELL.EQ.1) THEN
               COST = GPHAK(BETA)
            ELSEIF(ISHELL.EQ.2) THEN
               COST = GPHAL1(BETA)
            ELSEIF(ISHELL.EQ.3) THEN
               COST = GPHAL2(BETA)
            ELSEIF(ISHELL.EQ.4) THEN
               COST = GPHAL3(BETA)
            ENDIF
            PHI = TWOPI*RN03
            COSPHI = COS(PHI)
            SINPHI = SIN(PHI)
            SINT = SQRT((1.-COST)*(1.+COST))
            GKIN(1,NGKINE) = PEL*SINT*COSPHI
            GKIN(2,NGKINE) = PEL*SINT*SINPHI
            GKIN(3,NGKINE) = PEL*COST
            GKIN(4,NGKINE) = EEL
            GKIN(5,NGKINE) = 3.
            TOFD(NGKINE) = 0
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
C
C             ROTATE ELECTRON AND SCATTERED PHOTON INTO GEANT SYSTEM
C
            CALL GVROT(VECT(4),GKIN)
         ENDIF
         IF(ISHELL.NE.0) THEN
            ISHELL = ABS(ISHELL)
            IF(ZINT.GE.5.AND.POT(ISHELL).GT.MIN(CUTGAM,CUTELE)) THEN
C Generate shell decay mode
               IF(RN04.LE.Q(JFN+1+NSHELL+ISHELL)) THEN
                  IF(POT(ISHELL).LE.CUTGAM) THEN
                     DESTEP = DESTEP+POT(ISHELL)
                  ELSE
C Radiative shell decay
                     JS = JFN+1+2*NSHELL+ISHELL
                     JS = JPHFN+Q(JS)
                     NPOINT = Q(JS)
                     DO 40 I = 1,NPOINT
                        IF(RN05.LT.Q(JS+I)) THEN
                           TSEC = Q(JS+NPOINT+I)
                           IF(TSEC.GT.CUTGAM) THEN
                              NGKINE = NGKINE+1
                              PHI = TWOPI*RN06
                              COST = 2.*RN07-1.
                              COSPHI = COS(PHI)
                              SINPHI = SIN(PHI)
                              SINT = SQRT((1.-COST)*(1.+COST))
                              GKIN(1,NGKINE) = TSEC*SINT*COSPHI
                              GKIN(2,NGKINE) = TSEC*SINT*SINPHI
                              GKIN(3,NGKINE) = TSEC*COST
                              GKIN(4,NGKINE) = TSEC
                              GKIN(5,NGKINE) = 1.
                              TOFD(NGKINE) = 0.
                              GPOS(1,NGKINE) = VECT(1)
                              GPOS(2,NGKINE) = VECT(2)
                              GPOS(3,NGKINE) = VECT(3)
                           ELSE
                              DESTEP = DESTEP+ABS(TSEC)
                           ENDIF
C The following particle forces the energy conservation
                           TSEC = POT(ISHELL)-ABS(TSEC)
                           IF(TSEC.GT.CUTGAM) THEN
                              NGKINE = NGKINE+1
                              PHI = TWOPI*RN08
                              COST = 2.*RN09-1.
                              COSPHI = COS(PHI)
                              SINPHI = SIN(PHI)
                              SINT = SQRT((1.-COST)*(1.+COST))
                              GKIN(1,NGKINE) = TSEC*SINT*COSPHI
                              GKIN(2,NGKINE) = TSEC*SINT*SINPHI
                              GKIN(3,NGKINE) = TSEC*COST
                              GKIN(4,NGKINE) = TSEC
                              GKIN(5,NGKINE) = 1.
                              TOFD(NGKINE) = 0.
                              GPOS(1,NGKINE) = VECT(1)
                              GPOS(2,NGKINE) = VECT(2)
                              GPOS(3,NGKINE) = VECT(3)
                           ELSE
                              DESTEP = DESTEP+TSEC
                           ENDIF
                           GO TO 50
                        ENDIF
   40                CONTINUE
   50                CONTINUE
                  ENDIF
               ELSE
                  IF(POT(ISHELL).LE.CUTELE) THEN
                     DESTEP = DESTEP+POT(ISHELL)
                  ELSE
c Nonradiative decay
                     JS = JFN+1+3*NSHELL+ISHELL
                     JS = JPHFN+Q(JS)
                     NPOINT = Q(JS)
                     DO 60 I = 1,NPOINT
                        IF(RN05.LT.Q(JS+I)) THEN
                           TSEC = Q(JS+NPOINT+I)
                           IF(TSEC.GT.CUTELE) THEN
                              EEL=TSEC + EMASS
                              PEL=SQRT((TSEC+2.*EMASS)*TSEC)
                              NGKINE = NGKINE+1
                              PHI = TWOPI*RN06
                              COST = 2.*RN07-1.
                              COSPHI = COS(PHI)
                              SINPHI = SIN(PHI)
                              SINT = SQRT((1.-COST)*(1.+COST))
                              GKIN(1,NGKINE) = PEL*SINT*COSPHI
                              GKIN(2,NGKINE) = PEL*SINT*SINPHI
                              GKIN(3,NGKINE) = PEL*COST
                              GKIN(4,NGKINE) = EEL
                              GKIN(5,NGKINE) = 3.
                              TOFD(NGKINE) = 0.
                              GPOS(1,NGKINE) = VECT(1)
                              GPOS(2,NGKINE) = VECT(2)
                              GPOS(3,NGKINE) = VECT(3)
                           ELSE
                              DESTEP = DESTEP+ABS(TSEC)
                           ENDIF
C The following particle forces the energy conservation
                           TSEC = POT(ISHELL)-ABS(TSEC)
                           IF(TSEC.GT.CUTGAM) THEN
                              NGKINE = NGKINE+1
                              PHI = TWOPI*RN08
                              COST = 2.*RN09-1.
                              COSPHI = COS(PHI)
                              SINPHI = SIN(PHI)
                              SINT = SQRT((1.-COST)*(1.+COST))
                              GKIN(1,NGKINE) = TSEC*SINT*COSPHI
                              GKIN(2,NGKINE) = TSEC*SINT*SINPHI
                              GKIN(3,NGKINE) = TSEC*COST
                              GKIN(4,NGKINE) = TSEC
                              GKIN(5,NGKINE) = 1.
                              TOFD(NGKINE) = 0.
                              GPOS(1,NGKINE) = VECT(1)
                              GPOS(2,NGKINE) = VECT(2)
                              GPOS(3,NGKINE) = VECT(3)
                           ELSE
                              DESTEP = DESTEP+TSEC
                           ENDIF
                           GO TO 70
                        ENDIF
   60                CONTINUE
   70                CONTINUE
                  ENDIF
               ENDIF
            ELSE
               DESTEP = DESTEP+POT(ISHELL)
            ENDIF
         ENDIF
      ENDIF
      END

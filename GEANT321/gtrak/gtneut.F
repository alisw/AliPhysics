*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.24  by  S.Giani
*-- Author :
      SUBROUTINE GTNEUT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Neutral hadron type track. Computes step size and propagates *
C.    *    particle through step.                                      *
C.    *                                                                *
C.    *   ==>Called by : GTRACK                                        *
C.    *      Authors     R.Brun, F.Bruyant  *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#if defined(CERNLIB_DEBUG)
#include "geant321/gcunit.inc"
#endif
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
C.
C.    ------------------------------------------------------------------
*
* *** Particle below energy threshold ? Special treatment
*
      IF (GEKIN.LE.CUTNEU) THEN
         GEKIN  = 0.
         GETOT  = AMASS
         VECT(7)= 0.
         ISTOP  = 2
         NMEC = NMEC + 1
         LMEC(NMEC) = 30
         IF (IHADR.EQ.0.AND.AMASS.GT.0.) THEN
            IF (IDCAY.NE.0) THEN
               GAMMA  = GETOT/AMASS
               TOFG = TOFG +GAMMA*SUMLIF/CLIGHT
               SUMLIF = 0.
               IF (TOFG.GE.TOFMAX) THEN
                  ISTOP = 4
                  NMEC = NMEC + 1
                  LMEC(NMEC) = 22
                  NGKINE = 0
               ELSE
                  NMEC = NMEC + 1
                  LMEC(NMEC) = 5
                  ISTOP =1
                  CALL GDECAY
               ENDIF
            ENDIF
            GO TO 999
         ENDIF
         IPROC = 12
         GO TO 40
      ENDIF
*
* *** Compute step size
*
      IPROC  = 103
      STEP   = STEMAX
*
*  **   Step limitation due to hadron interaction ?
*
      IF (IHADR.GT.0) THEN
#if !defined(CERNLIB_USRJMP)
         CALL GUPHAD
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUPHAD)
#endif
         IF (SHADR.LT.STEP) THEN
            IF (SHADR.LE.0.) SHADR = PREC
            STEP  = SHADR
            IPROC = 12
         ENDIF
      ENDIF
*
*  **   Step limitation due to decay ?
*
      IF ((IDCAY.NE.0).AND.(AMASS.GT.0.)) THEN
         SDCAY = SUMLIF*VECT(7)/AMASS
         IF (SDCAY.LT.STEP) THEN
            STEP  = SDCAY
            IPROC = 5
         ENDIF
      ENDIF
*
      IF (STEP.LT.0.) STEP = 0.
*
*  **   Step limitation due to geometry ?
*
      IF (STEP.GE.SAFETY) THEN
         CALL GTNEXT
         IF (IGNEXT.NE.0) THEN
            STEP = SNEXT + PREC
            IPROC = 0
            INWVOL = 2
            NMEC = NMEC + 1
            LMEC(NMEC)  = 1
         ENDIF
*
*        Update SAFETY in stack companions, if any
         IF (IQ(JSTAK+3).NE.0) THEN
            DO 10 IST = IQ(JSTAK+3),IQ(JSTAK+1)
               JST = JSTAK +3 +(IST-1)*NWSTAK
               Q(JST+11) = SAFETY
   10       CONTINUE
            IQ(JSTAK+3) = 0
         ENDIF
*
      ELSE
         IQ(JSTAK+3) = 0
      ENDIF
*
* *** Linear transport
*
      IF (INWVOL.EQ.2) THEN
         DO 20 I = 1,3
            VECTMP  = VECT(I) +STEP*VECT(I+3)
            IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
               IF(VECT(I+3).NE.0.) THEN
                  VECTMP = VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*
     +            EPSMAC
                  IF(NMEC.GT.0) THEN
                     IF(LMEC(NMEC).EQ.104) NMEC=NMEC-1
                  ENDIF
                  NMEC=NMEC+1
                  LMEC(NMEC)=104
#if defined(CERNLIB_DEBUG)
                  WRITE(CHMAIL, 10000)
                  CALL GMAIL(0,0)
                  WRITE(CHMAIL, 10100) GEKIN, NUMED, STEP, SNEXT
                  CALL GMAIL(0,0)
10000 FORMAT(' Boundary correction in GTNEUT: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
               ENDIF
            ENDIF
            VECT(I) = VECTMP
   20    CONTINUE
      ELSE
         DO 30 I = 1,3
            VECT(I)  = VECT(I) +STEP*VECT(I+3)
   30    CONTINUE
      ENDIF
*
      SLENG = SLENG +STEP
*
* *** Update time of flight
*
      SUMLIF = SUMLIF -STEP*AMASS/VECT(7)
      TOFG = TOFG +STEP*GETOT/(VECT(7)*CLIGHT)
      IF (TOFG.GE.TOFMAX) THEN
         ISTOP = 4
         NMEC  = NMEC +1
         LMEC(NMEC) = 22
         GO TO 999
      ENDIF
*
* *** Update interaction probabilities
*
      IF (IHADR.GT.0) ZINTHA = ZINTHA*(1.-STEP/SHADR)
*
* *** apply the selected process if any
*
   40 IF (IPROC.EQ.0) GO TO 999
      NMEC = NMEC +1
      LMEC(NMEC) = IPROC
*
*  **   Hadron interaction ?
*
      IF (IPROC.EQ.12) THEN
#if !defined(CERNLIB_USRJMP)
         CALL GUHADR
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUHADR)
#endif
*          Check time cut-off for decays at rest
         IF (LMEC(NMEC).EQ.5) THEN
            TOFG   = TOFG +SUMLIF/CLIGHT
            SUMLIF = 0.
            IF (TOFG.GE.TOFMAX) THEN
               ISTOP = 4
               LMEC(NMEC) = 22
               NGKINE = 0
            ENDIF
         ENDIF
*
*  **   Decay ?
*
      ELSE IF (IPROC.EQ.5) THEN
         ISTOP = 1
         CALL GDECAY
      ENDIF
*                                                             END GTNEUT
  999 END

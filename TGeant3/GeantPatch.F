*CMZ :          02/02/99  18.09.13  by  Federico Carminati
*-- Author :    Federico Carminati   02/02/99
      SUBROUTINE GFLUCT(DEMEAN,DE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Subroutine to decide which method is used to simulate        *
C.    *   the straggling around the mean energy loss.                  *
C.    *                                                                *
C.    *                                                                *
C.    *   DNMIN:  <---------->1<-------->30<--------->50<--------->    *
C.    *                                                                *
C.    *   LOSS=2  :                                                    *
C.    *   STRA=0  <----------GLANDZ-------------------><--GLANDO-->    *
C.    *   LOSS=1,3:                                                    *
C.    *   STRA=0  <---------------------GLANDZ-------------------->    *
C.    *                                                                *
C.    *   STRA=1  <-----------PAI---------------------><--GLANDZ-->    *
C.    *                                                                *
C.    *   DNMIN :  an estimation of the number of collisions           *
C.    *            with energy close to the ionization energy          *
C.    *            (see PHYS333)                                       *
C.    *                                                                *
C.    *   Input  : DEMEAN (mean energy loss)                           *
C.    *   Output : DE   (energy loss in the current step)              *
C.    *                                                                *
C.    *    ==>Called by : GTELEC,GTMUON,GTHADR                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/pilot.h"
#undef CERNLIB_GEANT321_GCBANK_INC
#undef CERNLIB_GEANT321_GCLINK_INC
#include "geant321/gcbank.inc"
#undef CERNLIB_GEANT321_GCJLOC_INC
#include "geant321/gcjloc.inc"
#undef CERNLIB_GEANT321_GCONSP_INC
#include "geant321/gconsp.inc"
#undef CERNLIB_GEANT321_GCMATE_INC
#include "geant321/gcmate.inc"
#undef CERNLIB_GEANT321_GCCUTS_INC
#include "geant321/gccuts.inc"
#undef CERNLIB_GEANT321_GCKINE_INC
#include "geant321/gckine.inc"
#undef CERNLIB_GEANT321_GCMULO_INC
#include "geant321/gcmulo.inc"
#undef CERNLIB_GEANT321_GCPHYS_INC
#include "geant321/gcphys.inc"
#undef CERNLIB_GEANT321_GCTRAK_INC
#include "geant321/gctrak.inc"
*KEND.
      PARAMETER (EULER=0.577215,GAM1=EULER-1)
      PARAMETER (P1=.60715,P2=.67794,P3=.52382E-1,P4=.94753,
     +           P5=.74442,P6=1.1934)
      PARAMETER (DGEV=0.153536 E-3, DNLIM=50)
* These parameters are needed by M.Kowalski's fluctuation algorithm
      PARAMETER (FPOT=20.77E-9, EEND=10E-6, EEXPO=2.2)
      PARAMETER (XEXPO=-EEXPO+1, YEXPO=1/XEXPO)
* These parameters are needed by M.Kowalski's fluctuation algorithm
      DIMENSION RNDM(2)
      DE2(DPOT,RAN)=(DPOT**XEXPO*(1-RAN)+EEND**XEXPO*RAN)**YEXPO
      FLAND(X) = P1+P6*X+(P2+P3*X)*EXP(P4*X+P5)
      IF(STEP.LE.0) THEN
         DE=DEMEAN
      ELSE
         DEDX = DEMEAN/STEP
         POTI=Q(JPROB+9)
         IF(ISTRA.EQ.0.AND.(ILOSS.EQ.1.OR.ILOSS.EQ.3)) THEN
            CALL GLANDZ(Z,STEP,VECT(7),GETOT,DEDX,DE,POTI,Q(JPROB+10))
         ELSEIF (ILOSS.EQ.5) THEN
* This is Marek Kowalski's fluctuation algorithm, it works only when
* the step size has been limited to one ionisation on average
            CALL GRNDM(RNDM,1)
            DE=DE2(FPOT,RNDM(1))
*
         ELSE
* *** mean ionization potential (GeV)
*        POTI=16E-9*Z**0.9
            GAMMA = GETOT/AMASS
            BETA = VECT(7)/GETOT
            BET2 = BETA**2
* ***    low energy transfer
            XI = DGEV*CHARGE**2*STEP*DENS*Z/(A*BET2)
* ***    regime
* ***    ISTRA = 1 --> PAI + URBAN
            DNMIN = MIN(XI,DEMEAN)/POTI
            IF (ISTRA.EQ.0) THEN
               IF(DNMIN.GE.DNLIM) THEN
*  Energy straggling using Gaussian, Landau & Vavilov theories.
*  STEP   =  current step-length (cm)
*  DELAND =  DE/DX - <DE/DX>     (GeV)
*  Author      : G.N. Patrick
                  IF(STEP.LT.1.E-7)THEN
                     DELAND=0.
                  ELSE
*     Maximum energy transfer to atomic electron (GeV).
                     ETA = BETA*GAMMA
                     RATIO = EMASS/AMASS
*     Calculate Kappa significance ratio.
*                 EMAX=(2*EMASS*ETA**2)/(1+2*RATIO*GAMMA+RATIO**2)
*                 CAPPA = XI/EMAX
                     CAPPA = XI*(1+2*RATIO*GAMMA+RATIO**2)/(2*EMASS*
     +               ETA**2)
                     IF (CAPPA.GE.10.) THEN
*     +-----------------------------------+
*     I Sample from Gaussian distribution I
*     +-----------------------------------+
                        SIGMA = XI*SQRT((1.-0.5*BET2)/CAPPA)
                        CALL GRNDM(RNDM,2)
                        F1 = -2.*LOG(RNDM(1))
                        DELAND = SIGMA*SQRT(F1)*COS(TWOPI*RNDM(2))
                     ELSE
                        XMEAN = -BET2-LOG(CAPPA)+GAM1
                        IF (CAPPA.LT.0.01) THEN
                           XLAMX = FLAND(XMEAN)
*     +---------------------------------------------------------------+
*     I Sample lambda variable from Kolbig/Schorr Landau distribution I
*     +---------------------------------------------------------------+
*  10                   CALL GRNDM(RNDM,1)
*                       IF( RNDM(1) .GT. 0.980 ) GO TO 10
*                       XLAMB = GLANDR(RNDM(1))
*     +---------------------------------------------------------------+
*     I Sample lambda variable from James/Hancock Landau distribution I
*     +---------------------------------------------------------------+
   10                      CALL GLANDG(XLAMB)
                           IF(XLAMB.GT.XLAMX) GO TO 10
                        ELSE
*            +---------------------------------------------------------+
*            I Sample lambda variable (Landau not Vavilov) from        I
*            I Rotondi&Montagna&Kolbig Vavilov distribution            I
*            +---------------------------------------------------------+
                           CALL GRNDM(RNDM,1)
                           XLAMB = GVAVIV(CAPPA,BET2,RNDM(1))
                        ENDIF
*     Calculate DE/DX - <DE/DX>
                        DELAND = XI*(XLAMB-XMEAN)
                     ENDIF
                  ENDIF
                  DE = DEMEAN + DELAND
               ELSE
                  CALL GLANDZ(Z,STEP,VECT(7),GETOT,DEDX,DE,POTI,
     +            Q(JPROB+ 10))
               ENDIF
            ELSE IF (ISTRA.LE.2) THEN
               IF(DNMIN.GE.DNLIM) THEN
                  CALL GLANDZ(Z,STEP,VECT(7),GETOT,DEDX,DE,POTI,
     +            Q(JPROB+ 10))
               ELSE
                  NMEC = NMEC+1
                  LMEC(NMEC)=109
                  DE = GSTREN(GAMMA,DCUTE,STEP)
* ***   Add brem losses to ionisation
                  IF(ITRTYP.EQ.2) THEN
                     JBASE = LQ(JMA-1)+2*NEK1+IEKBIN
                     DE = DE +(1.-GEKRAT)*Q(JBASE)+GEKRAT*Q(JBASE+1)
                  ELSEIF(ITRTYP.EQ.5) THEN
                     JBASE = LQ(JMA-2)+NEK1+IEKBIN
                     DE = DE +(1.-GEKRAT)*Q(JBASE)+GEKRAT*Q(JBASE+1)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      END
*CMZ :          16/02/99  19.24.38  by  Federico Carminati
*CMZ :  2.03/01 28/08/98  09.33.11  by  Federico Carminati
*CMZ :  2.01/00 18/06/98  17.34.13  by  Federico Carminati
*CMZ :  2.00/05 28/05/98  19.04.13  by  Federico Carminati
*-- Author :
      SUBROUTINE GTREVE
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GTREVE                                                *
C.    *                                                                *
C.    *   Controls tracking of all particles belonging to the current  *
C.    *    event.                                                      *
C.    *                                                                *
C.    *   Called by : GUTREV, called by GTRIG                          *
C.    *   Authors   : R.Brun, F.Bruyant                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/pilot.h"
#undef CERNLIB_GEANT321_GCBANK_INC
#undef CERNLIB_GEANT321_GCLINK_INC
#include "geant321/gcbank.inc"
#undef CERNLIB_GEANT321_GCFLAG_INC
#include "geant321/gcflag.inc"
#undef CERNLIB_GEANT321_GCKINE_INC
#include "geant321/gckine.inc"
#undef CERNLIB_GEANT321_GCNUM_INC
#include "geant321/gcnum.inc"
#undef CERNLIB_GEANT321_GCSTAK_INC
#include "geant321/gcstak.inc"
#undef CERNLIB_GEANT321_GCTMED_INC
#include "geant321/gctmed.inc"
#undef CERNLIB_GEANT321_GCTRAK_INC
#include "geant321/gctrak.inc"
#undef CERNLIB_GEANT321_GCUNIT_INC
#include "geant321/gcunit.inc"
#include "sckine.inc"
*KEND.
*
      REAL UBUF(2)
      EQUIVALENCE (UBUF(1),WS(1))
      LOGICAL   BTEST
      DIMENSION PMOM(3),VPOS(3)
C.
C.    ------------------------------------------------------------------
      NTMSTO = 0
      NSTMAX = 0
      NALIVE = 0
*         Kick start the creation of the vertex
      VPOS(1)=0
      VPOS(2)=0
      VPOS(3)=0
      PMOM(1)=0
      PMOM(2)=0
      PMOM(3)=0
      IPART=1
      CALL GSVERT(VPOS,0,0,UBUF,0,NVTX)
      CALL GSKINE(PMOM,IPART,NVTX,UBUF,0,NT)
*
      MTRACK=-999
 10   MTROLD=MTRACK
      CALL RXGTRAK(MTRACK,IPART,PMOM,E,VPOS,TTOF)
      IF(MTROLD.LT.0) THEN
         MPRIMA=MTRACK
      ENDIF
      IF(MTRACK.LE.MPRIMA) THEN
         IF(ISWIT(4).GT.0.AND.MTRACK.GT.0) THEN
            IF(ISWIT(4).EQ.1.OR.MOD(MTRACK,ISWIT(4)).EQ.0) THEN
               WRITE(CHMAIL,10200) MTRACK
               CALL GMAIL(0,0)
            ENDIF
         ENDIF
         IF(MTROLD.GT.0) THEN
C --- Output root hits tree only for each primary MTRACK
            CALL RXOUTH
         ENDIF
      ENDIF
      IF(MTRACK.LE.0) GOTO 999
* Set the vertex
      JV=LQ(JVERTX-1)
      Q(JV + 1) = VPOS(1)
      Q(JV + 2) = VPOS(2)
      Q(JV + 3) = VPOS(3)
      Q(JV + 4) = TTOF
      Q(JV + 5) = 0
      Q(JV + 6) = 0
* Set the track
      JK=LQ(JKINE-1)
      Q(JK + 1) = PMOM(1)
      Q(JK + 2) = PMOM(2)
      Q(JK + 3) = PMOM(3)
      Q(JK + 4) = E
      Q(JK + 5) = IPART
      Q(JK + 6) = 1
* Now transport
C      CALL GPVERT(0)
C      CALL GPKINE(0)
* Normal Gtreve code
      NV = NVERTX
      DO 40  IV = 1,NV
* ***   For each vertex in turn ..
         JV = LQ(JVERTX-IV)
         NT = Q(JV+7)
         IF (NT.LE.0) GO TO 40
         TOFG   = Q(JV+4)
         SAFETY = 0.
         IF (NJTMAX.GT.0) THEN
            CALL GMEDIA (Q(JV+1), NUMED)
            IF (NUMED.EQ.0) THEN
               WRITE (CHMAIL, 10000) (Q(JV+I), I=1,3)
               CALL GMAIL (0, 0)
               GO TO 40
            ENDIF
            CALL GLSKLT
         ENDIF
*  **   Loop over tracks attached to current vertex
         DO 20  IT = 1,NT
            JV   = LQ(JVERTX-IV)
            ITRA = Q(JV+7+IT)
            IF (BTEST(IQ(LQ(JKINE-ITRA)),0)) GO TO 20
            CALL GFKINE (ITRA, VERT, PVERT, IPART, IVERT, UBUF, NWBUF)
            IF (IVERT.NE.IV) THEN
               WRITE (CHMAIL, 10100) IV, IVERT
               CALL GMAIL (0, 0)
               GO TO 999
            ENDIF
*   *      Store current track parameters in stack JSTAK
            CALL GSSTAK (2)
   20    CONTINUE
*  **   Start tracking phase
   30    IF (NALIVE.NE.0) THEN
            NALIVE = NALIVE -1
*   *      Pick-up next track in stack JSTAK, if any
*   *         Initialize tracking parameters
            CALL GLTRAC
            IF (NUMED.EQ.0) GO TO 30
*   *       Resume tracking
            CALL GUTRAK
            IF (IEOTRI.NE.0) GO TO 999
            GO TO 30
         ENDIF
*
   40 CONTINUE
      GOTO 10
*
10000 FORMAT (' GTREVE : Vertex outside setup, XYZ=',3G12.4)
10100 FORMAT (' GTREVE : Abnormal track/vertex connection',2I8)
10200 FORMAT (' GTREVE : Transporting primary track No ',I10)
*                                                             END GTREVE
  999 END

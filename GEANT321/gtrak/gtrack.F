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
*FCA :          17/05/99  16:21:12  by  Federico Carminati
*               Added the modifications of P.Nevski in MANY volumes 
*               force update of alternative list of many candidates
*-- Author :
      SUBROUTINE GTRACK
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Controls tracking of current particle,                   *
C.    *        up to end of track for sequential tracking mode, or     *
C.    *        through current volume for parallel tracking mode.      *
C.    *                                                                *
C.    *    ==>Called by : GUTRAK                                       *
C.    *       Authors   : R.Brun, F.Bruyant                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcparm.inc"
#include "geant321/gcsets.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcnum.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      COMMON/GCCHAN/LSAMVL
      LOGICAL LSAMVL
*
      DIMENSION CUTS(10),MECA(5,13)
      EQUIVALENCE (CUTS(1),CUTGAM),(MECA(1,1),IPAIR)
      SAVE PRECOR
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
C.
C.    ------------------------------------------------------------------
      ISTOP = 0
      EPSCUR = EPSMAC
      NSTOUT = 0
      INWOLD = 0
      LSAMVL = .FALSE.
*
* *** Check validity of tracking medium and material parameters
*
   10 IF (NUMED.NE.NUMOLD) THEN
         NUMOLD = NUMED
         IUPD   = 0
         JTM    = LQ(JTMED- NUMED)
         DO 20 I = 1,5
            NATMED(I) = IQ(JTM+I)
   20    CONTINUE
         NMAT     = Q(JTM + 6)
         ISVOL    = Q(JTM + 7)
         IFIELD   = Q(JTM + 8)
         FIELDM   = Q(JTM + 9)
         TMAXFD   = Q(JTM + 10)
         STEMAX   = Q(JTM + 11)
         DEEMAX   = Q(JTM + 12)
         EPSIL    = Q(JTM + 13)
         STMIN    = Q(JTM + 14)
         PRECOR   = MIN(0.1*EPSIL, 0.0010)
         IF (LQ(JTM).EQ.0) THEN
            IF (ISTPAR.NE.0) THEN
               DO 30 I = 1,10
                  CUTS(I) = Q(JTMED+I)
   30          CONTINUE
               DO 40 I = 1,13
                  MECA(1,I) = Q(JTMED+10+I)
   40          CONTINUE
               ILABS = Q(JTMED+10+21)
               ISYNC = Q(JTMED+10+22)
               ISTRA = Q(JTMED+10+23)
               ISTPAR = 0
            ENDIF
         ELSE
            JTMN = LQ(JTM)
            DO 50 I = 1,10
               CUTS(I) = Q(JTMN+I)
   50       CONTINUE
            DO 60 I = 1,13
               MECA(1,I) = Q(JTMN+10+I)
   60       CONTINUE
            ILABS = Q(JTMN+10+21)
            ISYNC = Q(JTMN+10+22)
            ISTRA = Q(JTMN+10+23)
            ISTPAR = 1
         ENDIF
*
         JMA   = LQ(JMATE-NMAT)
         JPROB = LQ(JMA-4)
         JMIXT = LQ(JMA-5)
         DO 70 I = 1,5
            NAMATE(I) = IQ(JMA+I)
   70    CONTINUE
         A    = Q(JMA +6)
         Z    = Q(JMA +7)
         DENS = Q(JMA +8)
         RADL = Q(JMA +9)
         ABSL = Q(JMA +10)
         IF(IQ(JTM-2).GE.3.AND.LQ(JTM-3).NE.0.AND.ITCKOV.NE.0.AND.
     +      LQ(LQ(JTM-3)-3).NE.0.AND.Z.GE.1.) THEN
*
* ***  In this tracking medium Cerenkov photons are generated and
* ***  tracked. Set to 1 the corresponding flag.
*
            IMCKOV = 1
         ELSE
            IMCKOV = 0
         ENDIF
*
*
*  **   Update precomputed quantities
*
         IMULL = IMULS
         IF (ILOSS.LE.0) THEN
            DEEMAX = 0.
            ILOSL = 0
         ELSEIF (DEEMAX.GT.0.) THEN
            ILOSL = ILOSS
         ELSE
            ILOSL = 0
         ENDIF
      ENDIF
*
      IF(LSAMVL) THEN
*
*       If now the particle is entering in the same volume where
*       it was exiting from last step, and if it has done this for
*       more than 5 times, we decrease the precision of tracking
         NSTOUT=NSTOUT+1
         IF(MOD(NSTOUT,5).EQ.0) THEN
            EPSCUR=NSTOUT*EPSMAC
*            WRITE(CHMAIL,10000)ITRA,ISTAK,NTMULT,NAPART
*10000          FORMAT(' *** GTRACK *** Boundary loop: track ',
*     +         I4,' stack ',I4,' NTMULT ',I5,1X,5A4)
*            CALL GMAIL(1,0)
*            WRITE (CHMAIL,10250) IEVENT,IDEVT,(NRNDM(I),I = 1,2)
*            CALL GMAIL(0,0)
*            WRITE(CHMAIL,10100) EPSCUR
*10100          FORMAT('                Precision now set to ',G10.3)
*            CALL GMAIL(0,1)
         ENDIF
      ELSE
         NSTOUT = 0
         EPSCUR = EPSMAC
      ENDIF
*
      INWVOL = 1
*
* *** Compute SET and DET number if volume is sensitive
*
      IF (JSET.GT.0) THEN
         IF(ISVOL.GT.0) THEN
            CALL GFINDS
         ELSE
            IHSET = 0
            IHDET = 0
            ISET = 0
            IDET = 0
            IDTYPE = 0
            NVNAME = 0
         ENDIF
      ENDIF
*
*    Clear step dependent variables
*
   80 NMEC   = 0
      STEP   = 0.
      DESTEL = 0.
      DESTEP = 0.
      NGKINE = 0
      NGPHOT = 0.
      IGNEXT = 0
      INWOLD = INWVOL
      PREC   = MAX(PRECOR,MAX(ABS(VECT(1)),ABS(VECT(2)),
     +                        ABS(VECT(3)),SLENG)*EPSCUR)
*
*     Give control to user at entrance of volume (INWVOL=1)
*
      IF (INWVOL.EQ.1) THEN
#if !defined(CERNLIB_USRJMP)
         CALL GUSTEP
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUSTEP)
#endif
         IF (ISTOP.NE.0) GO TO 999
         INWVOL = 0
      ENDIF
*
* *** Propagate particle up to next volume boundary or end of track
*
      INGOTO = 0
      NLEVIN = NLEVEL
      IF (IPARAM.NE.0) THEN
         IF (GEKIN.LE.PACUTS(ITRTYP)) THEN
            NMEC = NMEC+1
            LMEC(NMEC) = 26
            ISTOP = 2
#if !defined(CERNLIB_USRJMP)
            CALL GUPARA
#endif
#if defined(CERNLIB_USRJMP)
            CALL JUMPT0(JUPARA)
#endif
            GO TO 90
         ENDIF
      ENDIF
      IF      (ITRTYP.EQ.1) THEN
         CALL GTGAMA
      ELSE IF (ITRTYP.EQ.2) THEN
         CALL GTELEC
      ELSE IF (ITRTYP.EQ.3) THEN
         CALL GTNEUT
      ELSE IF (ITRTYP.EQ.4) THEN
         CALL GTHADR
      ELSE IF (ITRTYP.EQ.5) THEN
         CALL GTMUON
      ELSE IF (ITRTYP.EQ.6) THEN
         CALL GTNINO
      ELSE IF (ITRTYP.EQ.7) THEN
         CALL GTCKOV
      ELSE IF (ITRTYP.EQ.8) THEN
         CALL GTHION
      ENDIF
      IF(JGSTAT.NE.0) CALL GFSTAT(10+ITRTYP)
      STLOSS=STEP
*
*     Check for possible endless loop
*
   90 NSTEP = NSTEP +1
      IF (NSTEP.GT.MAXNST) THEN
         IF (ISTOP.EQ.0) THEN
            ISTOP = 99
            NMEC  = NMEC +1
            LMEC(NMEC) = 30
            WRITE(CHMAIL,10200) MAXNST
            CALL GMAIL(1,0)
            WRITE(CHMAIL,10250) IEVENT,IDEVT,(NRNDM(I),I = 1,2)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10300)ITRA,ISTAK,NTMULT,(NAPART(I),I=1,5),
     +      TOFG*1.E9
            CALL GMAIL(0,1)
10200       FORMAT(' *** GTRACK *** More than ',I6,
     +             ' steps, tracking abandoned!')
10250       FORMAT('                IEVENT=',I7,' IDEVT=',I7,
     +             ' Random Seeds = ',I10,2X,I10)
10300       FORMAT('                Track ',I5,' stack ',I5,' NTMULT ',
     +             I5,1X,5A4,1X,'Time of flight ',F10.3,' ns')
         ENDIF
      ENDIF
*
* *** Give control to user at end of each tracking step
*
      SAFETY = SAFETY -STEP
#if !defined(CERNLIB_USRJMP)
      CALL GUSTEP
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUSTEP)
#endif
*
      IF (ISTOP.NE.0) GO TO 999
*
*      Renormalize direction cosines
*
      PMOM = SQRT(VECT(4)**2+VECT(5)**2+VECT(6)**2)
      IF(PMOM.GT.0.) THEN
         CMOD = 1./PMOM
         VECT(4) = VECT(4)*CMOD
         VECT(5) = VECT(5)*CMOD
         VECT(6) = VECT(6)*CMOD
      ENDIF
*                                        force update of alternatives:  
      IF (INWVOL.EQ.0) then
         if (GONLY(NLEVEL).eq.0.and.Safety.le.0.and.Iswit(10).ge.0)
     +        CALL GTMEDI(VECT,NMED)
         GO TO 80
      endif
*
      IF (NJTMAX.GT.0) THEN
         CALL GSTRAC
         IF (NLEVIN.EQ.0) GO TO 100
         GO TO 999
      ELSE
         IF (NLEVIN.GE.NLEVEL) THEN
            INFROM = 0
         ELSE
            IF (NLEVIN.EQ.0) GO TO 100
            INFROM = LINDEX(NLEVIN+1)
         ENDIF
         IF (NLEVIN.NE.NLEVEL) INGOTO = 0
         NLEVEL = NLEVIN
*
         CALL GTMEDI (VECT, NUMED)
         IF (NUMED.NE.0) THEN
            SAFETY = 0.
            GO TO 10
         ENDIF
      ENDIF
*
*     Track outside setup, give control to user (INWVOL=3)
*
  100 INWVOL = 3
      ISTOP  = 1
      ISET   = 0
      IDET   = 0
      NMEC   = 0
      STEP   = 0.
      DESTEL = 0.
      DESTEP = 0.
      NGKINE = 0
      NLCUR  = NLEVEL
      NLEVEL = 1
#if !defined(CERNLIB_USRJMP)
      CALL GUSTEP
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUSTEP)
#endif
      NLEVEL = NLCUR
*                                                             END GTRACK
  999 END

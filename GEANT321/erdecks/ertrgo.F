*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERTRGO
*
C.    ******************************************************************
C.    *                                                                *
C.    *       Perform the tracking of the track                        *
C.    *       Track parameters are in VECT                             *
C.    *                                                                *
C.    *    ==>Called by : ERTRAK                                       *
C.    *       Original routines : GTRACK + GTVOL                       *
C.    *       Authors   M.Maire, E.Nagy  *********                     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gccuts.inc"
#include "geant321/gconst.inc"
#include "geant321/gcphys.inc"
#include "geant321/gckine.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcmate.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
#include "geant321/ertrio.inc"
#include "geant321/erwork.inc"
      COMMON/GCCHAN/LSAMVL
      LOGICAL LSAMVL
*
*
      DIMENSION CUTS(10),MECA(5,12)
      EQUIVALENCE (CUTS(1),CUTGAM),(MECA(1,1),IPAIR)
      DIMENSION NAMIN(15),NUMIN(15),NAMOUT(15),NUMOUT(15)
*
      SAVE PRECOR,NSTOUT
#if (!defined(CERNLIB_SINGLE))&&(!defined(CERNLIB_IBM))
      PARAMETER (EPSMAC=5.E-6)
#endif
#if (!defined(CERNLIB_SINGLE))&&(defined(CERNLIB_IBM))
      PARAMETER (EPSMAC=5.E-5)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
C.
C.    ------------------------------------------------------------------
*
      NSTOUT = 0
      EPSCUR = EPSMAC
      LSAMVL = .FALSE.
      SLENG  = 0.
      ISTOP  = 0
      NUMED  = 0
      NUMOLD = 0
      IUPD   = 0
      NMEC   = 0
      INGOTO = 0
      INFROM = 0
      SAFETY = 0.
      MXNSTP = 1000
      NSTEP  = 0
*
      CALL GMEDIA(VECT,NUMED)
      IF (NUMED.EQ.0) GO TO 200
*
* *** Come back here each time we enter into a new volume
*
   10 CONTINUE
*
* *** Get tracking medium and material parameters
      IF (NUMED.NE.NUMOLD) THEN
         NUMOLD = NUMED
         IUPD   = 0
         JTM = LQ(JTMED- NUMED)
         DO 20 I=1,5
            NATMED(I)=IQ(JTM+I)
  20     CONTINUE
         NMAT   = Q(JTM + 6)
         ISVOL  = Q(JTM + 7)
         IFIELD = Q(JTM + 8)
         FIELDM = Q(JTM + 9)
         TMAXFD = Q(JTM + 10)
         DMAXMS = Q(JTM + 11)
         DEEMAX = Q(JTM + 12)
         EPSIL  = Q(JTM + 13)
         STMIN  = Q(JTM + 14)
         PRECOR   = MIN(0.1*EPSIL, 0.0010)
*
         IF(LQ(JTM).EQ.0)THEN
            IF(ISTPAR.NE.0)THEN
               DO 30 I=1,10
                  CUTS(I)=Q(JTMED+I)
  30           CONTINUE
               DO 40 I=1,12
                  MECA(1,I)=Q(JTMED+10+I)
  40           CONTINUE
               ISTPAR=0
            ENDIF
         ELSE
            JTMN=LQ(JTM)
            DO 50 I=1,10
               CUTS(I)=Q(JTMN+I)
  50        CONTINUE
            DO 60 I=1,12
               MECA(1,I)=Q(JTMN+10+I)
  60        CONTINUE
            ILABS = Q(JTMN+10+21)
            ISYNC = Q(JTMN+10+22)
            ISTRA = Q(JTMN+10+23)
            ISTPAR=1
         ENDIF
*
         JMA  = LQ(JMATE- NMAT)
         JPROB=LQ(JMA-4)
         JMIXT=LQ(JMA-5)
         DO 70 I=1,5
  70     NAMATE(I)=IQ(JMA+I)
         A    = Q(JMA + 6)
         Z    = Q(JMA + 7)
         DENS = Q(JMA + 8)
         RADL = Q(JMA + 9)
         ABSL = Q(JMA + 10)
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
            WRITE(CHMAIL,10000)ITRA,ISTAK,NTMULT,NAPART
10000          FORMAT(' *** ERTRGO *** Boundary loop: track ',
     +         I4,' stack ',I4,' NTMULT ',I5,1X,5A4)
            CALL GMAIL(1,0)
            WRITE(CHMAIL,10100) EPSCUR
10100          FORMAT('                Precision now set to ',G10.3)
            CALL GMAIL(0,1)
         ENDIF
      ELSE
         NSTOUT = 0
         EPSCUR = EPSMAC
      ENDIF
*
* *** Initialize magnetic field for EMC package
      HI(1) = 0.
      HI(2) = 0.
      HI(3) = 0.
      IF (IFIELD.EQ.3) THEN
         HI(3) = FIELDM
      ELSEIF (IFIELD.NE.0) THEN
         CALL GUFLD (VECT, HI)
      ENDIF
*
* *** Control given to user at entrance of volume (INWVOL=1)
      INWVOL = 1
      NMEC   = 1
      LMEC(1) = 29
      STEP   = 0.
      DESTEP = 0.
*
      IF((LEVOLU).AND.(SLENG.GT.0.)) THEN
         IMEC = 0
         CALL EVOLIO(NVLIN,NAMIN,NUMIN,NVLOUT,NAMOUT,NUMOUT)
         DO 80 IPR =1,NEPRED
            NAMPR = NAMEER(IPR)
            NUMPR = NUMVER(IPR)
            IOVPR = IOVLER(IPR)
            IF (IOVPR.EQ.1) THEN
               IV = IUCOMP(NAMPR,NAMIN ,NVLIN )
               IF (IV.NE.0) THEN
                  IF (NUMPR.EQ.0) NUMPR = NUMIN (IV)
                  IF (NUMPR.EQ.NUMIN (IV)) THEN
                     NMEC = NMEC + 1
                     LMEC(NMEC) = 27
                     INLIST = IPR
                     CALL ERSTOR
                  ENDIF
               ENDIF
            ELSE IF (IOVPR.EQ.2) THEN
               IV = IUCOMP(NAMPR,NAMOUT,NVLOUT)
               IF (IV.NE.0) THEN
                  IF (NUMPR.EQ.0) NUMPR = NUMOUT(IV)
                  IF (NUMPR.EQ.NUMOUT(IV)) THEN
                     NMEC = NMEC + 1
                     LMEC(NMEC) = 27
                     INLIST = IPR
                     CALL ERSTOR
                  ENDIF
               ENDIF
            ENDIF
*
   80    CONTINUE
*
      ENDIF
*
      CALL EUSTEP
      IF (ISTOP.NE.0) GO TO 999
*
* *** Particle is propagated up to the next volume boundary
*
      INWVOL=0
*
* *** Come back here after each step in the same volume
 100  IGNEXT = 0
      INGOTO = 0
      NLEVIN = NLEVEL
      NMEC   = 0
      STEP   = 0.
      DESTEP = 0.
      DEDX2  = 0.
      PREC   = MAX(PRECOR,MAX(ABS(VECT(1)),ABS(VECT(2)),
     +                        ABS(VECT(3)),SLENG)*EPSCUR)
*
      IF(CHARGE.NE.0.) THEN
         CALL ERTRCH
      ELSE
         CALL ERTRNT
      ENDIF
*
      NSTEP  = NSTEP + 1
      IF (NSTEP.GT.MXNSTP) THEN
         ISTOP = 99
         NMEC  = NMEC + 1
         LMEC(NMEC) = 30
      ENDIF
*
      SAFETY = SAFETY - STEP
      TLRAD  = TLRAD  + STEP/RADL
      TLGCM2 = TLGCM2 + STEP*DENS
*
* *** Give control to user after each tracking step
      CALL EUSTEP
*
      IF(ISTOP.NE.0) GO TO 999
*
* *** Renormalize direction cosines
      CMOD = 1./SQRT(VECT(4)**2 + VECT(5)**2 + VECT(6)**2)
      VECT(4) = VECT(4)*CMOD
      VECT(5) = VECT(5)*CMOD
      VECT(6) = VECT(6)*CMOD
*
      IF (INWVOL.EQ.0) GO TO 100
*
* *** Particle is leaving the volume (INWVOL=2) :
*
*     Save the current volume's tree before leaving the volume
      IF(LEVOLU) CALL EVOLIO(NVLIN,NAMIN,NUMIN,NVLOUT,NAMOUT,NUMOUT)
*
*     find the new volume
      IF (NLEVIN.GE.NLEVEL) THEN
         INFROM = 0
      ELSE
         IF (NLEVIN.EQ.0) GO TO 200
         INFROM = LINDEX (NLEVIN+1)
      ENDIF
      IF (NLEVIN.NE.NLEVEL) INGOTO = 0
      NLEVEL = NLEVIN
*
      CALL GTMEDI (VECT,NUMED)
      IF (NUMED.NE.0) THEN
         SAFETY = 0.
         GO TO 10
      ENDIF
*
* *** Track outside setup, give control to user (INWVOL=3)
  200 INWVOL= 3
      ISTOP = 1
      NMEC  = NMEC + 1
      LMEC(NMEC) = 30
      CALL EUSTEP
 999  CONTINUE
      ILOSL = 0
*
      END

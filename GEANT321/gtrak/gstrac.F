*
* $Id$
*
* $Log$
* Revision 1.2  1996/02/27 10:30:56  ravndal
* Correct interaction length for heavy ions
*
* Revision 1.1.1.1  1995/10/24 10:21:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 13/12/94  15.23.45  by  S.Giani
*-- Author :
      SUBROUTINE GSTRAC
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       SUBR. GSTRAC                                             *
C.    *                                                                *
C.    *  Stores in stack JTRACK the information for current track      *
C.    *   segment at exit of current Volume/Medium.                    *
C.    *                                                                *
C.    *   Called by : GTRACK                                           *
C.    *   Authors   : S.Banerjee, F.Bruyant                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
      REAL      XC(3)
C.
C.    ------------------------------------------------------------------
*
* *** Find where in the tracking skeleton to enter the current track
*
      IF (NLEVIN.EQ.NLEVEL) THEN
*       (case when GINVOL has not been called)
         JVO = LQ(JVOLUM-LVOLUM(NLEVIN))
         IF (INGOTO.GT.0) THEN
*
*  **      Point is in content predicted by GTNEXT, go one level down
*
            NLEVIN = NLEVIN +1
            INFROM = 0
            JIN    = LQ(JVO-INGOTO)
            IVOT   = Q(JIN+2)
            LVOLUM(NLEVIN) = IVOT
            LINDEX(NLEVIN) = INGOTO
            LINMX (NLEVIN) = Q(JVO+3)
*
*  **      Prepare the translation and rotation matrices if necessary
*
            JSKL = LQ(JSKLT-NLEVIN)
            IF (NLEVIN.GT.2) THEN
               IOFF = IQ(JSKL-3)
               DO 29 ILEV = 1, NLEVEL
                  IF (IQ(JSKL+IOFF+ILEV).NE.LINDEX(ILEV)) GO TO 30
   29          CONTINUE
            ENDIF
            JSK  = LQ(JSKL-INGOTO)
            IF (IQ(JSK+1).GT.0) GO TO 100
   30       IROTT = Q(JIN+4)
            CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), Q(JIN+5),
     +                   IROTT, GTRAN(1,NLEVIN), GRMAT(1,NLEVIN))
            GO TO 100
         ELSE
*
*  **      otherwise, go one level up
*
            NLEVIN = NLEVIN -1
*
         ENDIF
*
      ELSE IF (NLEVIN.GT.NLEVEL) THEN
         INFROM = 0
         GO TO 100
*
      ELSE IF (NLEVIN.LT.0) THEN
*       (case when entering a dominant overlaping volume)
         NLEVIN = -NLEVIN
         INFROM = LINDEX(NLEVIN+1)
         GO TO 100
      ENDIF
*
*  **  Track has left current volume, check levels up
*
   80 IF (NLEVIN.EQ.0) GO TO 999
*
      IF (GRMAT(10,NLEVIN).EQ.0.) THEN
         DO 88 I = 1,3
            XC(I) = VECT(I) -GTRAN(I,NLEVIN)
   88    CONTINUE
      ELSE
C       (later, code in line)
         CALL GTRNSF (VECT, GTRAN(1,NLEVIN), GRMAT(1,NLEVIN), XC)
      ENDIF
*
      JVO  = LQ(JVOLUM-LVOLUM(NLEVIN))
      JPAR = LQ(JGPAR-NLEVIN)
      CALL GINME (XC, Q(JVO+2), Q(JPAR+1), IYES)
      IF (IYES.NE.0) THEN
         INFROM = LINDEX(NLEVIN+1)
      ELSE
         NLEVIN = NLEVIN -1
         GO TO 80
      ENDIF
*
* *** Allocate last 'garbaged' area if any, otherwise first 'free' one
*
  100 IF (NJGARB.NE.0) THEN
         NCUR   = NJGARB
         LCUR   = JTRACK +(NCUR-1)*NWTRAC
         NJGARB = IQ(LCUR+1)
      ELSE
         NCUR   = NJFREE
         LCUR   = JTRACK +(NCUR-1)*NWTRAC
         NJFREE = NCUR +1
      ENDIF
*
* *** Link allocated area to relevant chain in JSKLT structure
*
      JSKL = LQ(JSKLT-NLEVIN)
      IQ(LCUR+1) = IQ(JSKL+LINDEX(NLEVIN))
      IQ(JSKL+LINDEX(NLEVIN)) = NCUR
*
* *** Store information for current track segment in stack JTRACK
*
      IQ(LCUR+2) = 0
      IQ(LCUR+3) = NTMULT
      IQ(LCUR+4) = ITRA
      IQ(LCUR+5) = ISTAK
      IQ(LCUR+6) = IPART
      IQ(LCUR+7) = NSTEP
*free IQ(LCUR+8) = IDECAD
      IQ(LCUR+9) = IEKBIN
      IQ(LCUR+10)= ISTORY
      IQ(LCUR+11)= INFROM
*
      IPCUR = LCUR +NWINT
      DO 109 I = 1,7
         Q(IPCUR+I)  = VECT(I)
  109 CONTINUE
      Q(IPCUR+8)  = GEKIN
      Q(IPCUR+9)  = SLENG
      Q(IPCUR+10) = GEKRAT
      Q(IPCUR+11) = TOFG
      Q(IPCUR+12) = UPWGHT
*
      IPCUR = IPCUR +NWREAL
      IF (ITRTYP.EQ.1) THEN
*     Photons
         Q(IPCUR+1) = ZINTPA
         Q(IPCUR+2) = ZINTCO
         Q(IPCUR+3) = ZINTPH
         Q(IPCUR+4) = ZINTPF
         Q(IPCUR+5) = ZINTRA
      ELSE IF (ITRTYP.EQ.2) THEN
*     Electrons
         Q(IPCUR+1) = ZINTBR
         Q(IPCUR+2) = ZINTDR
         Q(IPCUR+3) = ZINTAN
      ELSE IF (ITRTYP.EQ.3) THEN
*     Neutral hadrons
         Q(IPCUR+1) = SUMLIF
         Q(IPCUR+2) = ZINTHA
      ELSE IF (ITRTYP.EQ.4) THEN
*     Charged hadrons
         Q(IPCUR+1) = SUMLIF
         Q(IPCUR+2) = ZINTHA
         Q(IPCUR+3) = ZINTDR
      ELSE IF (ITRTYP.EQ.5) THEN
*     Muons
         Q(IPCUR+1) = SUMLIF
         Q(IPCUR+2) = ZINTBR
         Q(IPCUR+3) = ZINTPA
         Q(IPCUR+4) = ZINTDR
         Q(IPCUR+5) = ZINTMU
      ELSE IF (ITRTYP.EQ.7) THEN
*     Cerenkov photons
         Q(IPCUR+1) = ZINTLA
      ELSE IF (ITRTYP.EQ.8) THEN
*     Ions
         Q(IPCUR+1) = ZINTHA
         Q(IPCUR+2) = ZINTDR
      ENDIF
*
* *** Take care of the skeleton
*
      IF (NLEVIN.GT.NLDOWN) THEN
         NLDOWN = NLEVIN
         JSKL   = LQ(JSKLT-NLDOWN)
*
*  **    Clear skeleton at lowest level if necessary
*
         JOFF   = JSKL + IQ(JSKL-3)
         DO 229 ILEV = 1, NLDOWN-1
            IF (IQ(JOFF+ILEV).EQ.LINDEX(ILEV)) GO TO 229
            NINSK = LINMX(NLDOWN)
            DO 209 IN = 1, NINSK
               JSK  = LQ(JSKL-IN)
               IQ(JSK+1) = 0
  209       CONTINUE
            DO 219 I = ILEV, NLDOWN-1
               IQ(JOFF+I) = LINDEX(I)
  219       CONTINUE
            GO TO 230
  229    CONTINUE
      ENDIF
*
*  ** Fill up the skeleton at NLDOWN
*
  230 IF (NLEVIN.GT.NLEVEL) THEN
         JSKL = LQ(JSKLT-NLDOWN)
         JSK  = LQ(JSKL-LINDEX(NLDOWN))
         IF (IQ(JSK+1).LE.0) THEN
            LQ(JSK-1) = LQ(JGPAR-NLDOWN)
            IQ(JSK+1) = IQ(JGPAR+NLDOWN)
            IQ(JSK+2) = LVOLUM(NLDOWN)
            DO 239 I = 1, 3
               Q(JSK+2+I) = GTRAN(I,NLDOWN)
  239       CONTINUE
            DO 249 I = 1, 10
               Q(JSK+5+I) = GRMAT(I,NLDOWN)
  249       CONTINUE
         ENDIF
*
      ENDIF
*
* *** Update NALIVE and test if tracking stack is full
*
      NALIVE = NALIVE + 1
      IF (NALIVE-IQ(JSTAK+1).GE.NJTMAX) THEN
         WRITE (CHMAIL, 1001)
         CALL GMAIL (0, 0)
         NJTMAX = -NJTMAX
         NLVSAV = NLEVEL
         DO 309 I = 2,NLDOWN
            LINSAV(I) = LINDEX(I)
            LMXSAV(I) = LINMX(I)
  309    CONTINUE
      ENDIF
*
 1001 FORMAT (' GSTRAC : Stack JTRACK full. Inhibit parallel tracking')
*                                                             END GSTRAC
  999 END
 

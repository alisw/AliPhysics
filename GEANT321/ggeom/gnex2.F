*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLD)
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GNEXT (X, SNEXT, SAFETY)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   SUBR. GNEXT (X, SNEXT, SAFETY)                               *
C.    *                                                                *
C.    *   Computes SNEXT and SAFETY                                    *
C.    *     SNEXT  (output) : distance to closest boundary             *
C.    *                      from point X(1-3) along X(4-6)            *
C.    *     SAFETY (output) : shortest distance to any boundary        *
C.    *                                                                *
C.    *   Called by : User                                             *
C.    *   Authors   : S.Banerjee, R.Brun, F.Bruyant                    *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcvolu.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      REAL    X(6), X0(3), XC(6), XT(6)
      INTEGER IDTYP(3,12)
      LOGICAL BTEST
      SAVE IDTYP
C.
      DATA  IDTYP / 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 2, 3, 1,
     +              2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 4, 3, 1, 1, 1,
     +              2, 3, 1, 2, 3, 1/
C.
C.    ------------------------------------------------------------------
*
* *** Transform current point and direction into local reference system
*
      IF (GRMAT(10,NLEVEL).EQ.0.) THEN
         DO 19 I = 1,3
            XC(I)   = X(I) -GTRAN(I,NLEVEL)
            XC(I+3) = X(I+3)
   19    CONTINUE
      ELSE
*       (later, code in line)
         CALL GTRNSF (X, GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), XC)
         CALL GROT (X(4), GRMAT(1,NLEVEL), XC(4))
      ENDIF
*
* *** Compute distance to boundaries
*
      SNEXT  = BIG
      SAFETY = BIG
      JVO    = LQ(JVOLUM-LVOLUM(NLEVEL))
      ISH    = Q(JVO+2)
      IF (Q(JVO+3).EQ.0.) GO TO 300
      NIN = Q(JVO+3)
      IF (NIN.LT.0) GO TO 200
*
* *** Case with contents positioned
*
      ISEARC = Q(JVO+1)
      IF (ISEARC.GE.-1) GO TO 120
*
*  ** Contents are ordered by (dynamic) GSORD, select neighbours
*
      JSB = LQ(LQ(JVO-NIN-1))
      IAX = Q(JSB+1)
      NSB = Q(JSB+2)
      IF (IAX.LE.3) THEN
         CX  = XC(IAX)
         INC = SIGN(1., XC(IAX+3))
      ELSE
         CALL GFCOOR (XC, IAX, CX)
         IF (IAX.LE.5) THEN
            DR = XC(1)*XC(4) +XC(2)*XC(5)
            IF (IAX.EQ.5) DR = DR +XC(3)*XC(6)
            INC = SIGN(1., DR)
         ELSE IF (IAX.EQ.6) THEN
            INC = SIGN(1., XC(1)*XC(5)-XC(2)*XC(4))
         ELSE
            INC = SIGN(1., XC(3)*(XC(1)*XC(4)+XC(2)*XC(5))
     +                    -XC(6)*(XC(1)*XC(1)+XC(2)*XC(2)))
         ENDIF
      ENDIF
      IDIV = LOCATF (Q(JSB+3), NSB, CX)
      IF (IDIV.LT.0) IDIV = -IDIV
      IF (IAX.NE.6) THEN
         IF (IDIV.EQ.0) THEN
            IF (INC.LT.0.AND.IAX.LE.3) THEN
               SAFETY = Q(JSB+3) -CX
               GO TO 300
            ENDIF
            IDIV = 1
         ELSE IF (IDIV.EQ.NSB) THEN
            IF (INC.GT.0.AND.IAX.NE.7) THEN
               SAFETY = CX -Q(JSB+2+NSB)
               GO TO 300
            ENDIF
            IDIV = NSB -1
         ELSE
            IF (IAX.NE.7) THEN
               IF (INC.GT.0) THEN
                  SAFETY = CX -Q(JSB+2+IDIV)
               ELSE
                  SAFETY = Q(JSB+3+IDIV) -CX
               ENDIF
            ELSE
               SAFETY = 0.
            ENDIF
         ENDIF
      ELSE IF (IAX.EQ.6) THEN
         IF (IDIV.EQ.0) IDIV = NSB
         SAFETY = 0.
      ENDIF
*
      IDIVL = 0
      IDIVB = 0
      JSC0  = LQ(JVO-NIN-2)
  110 NCONT = IQ(JSC0+IDIV)
*
*  ** Loop over (selected) contents
*
      IF (NCONT.EQ.0) THEN
         IF (IDIV.EQ.IDIVL) GO TO 400
         IDIV = IDIV +INC
*      (following statement for IAX=6, when division NSB is empty)
         IF (IDIV.GT.NSB) IDIV = 1
         GO TO 110
      ELSE
         ICONT = 1
         JSCV = LQ(JSC0-IDIV)
         GO TO 140
      ENDIF
*
  120 JNEAR = LQ(JVO-NIN-1)
      IF (ISEARC.GT.0) THEN
#if !defined(CERNLIB_USRJMP)
         CALL GUNEAR (ISEARC, 2, XC, JNEAR)
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT4(JUNEAR,ISEARC, 2, XC, JNEAR)
#endif
         IF (IQ(JNEAR+1).EQ.0) GO TO 300
      ENDIF
      JNEAR = JNEAR +1
      NNEAR = IQ(JNEAR)
      IF (IQ(JNEAR+1).NE.0) THEN
         INEAR = 1
      ELSE
         INEAR = 2
      ENDIF
*
  130 IN = IQ(JNEAR+INEAR)
      GO TO 150
*
  140 IN = IQ(JSCV+ICONT)
*
  150 IF(IN.LE.0)GO TO 300
      JIN   = LQ(JVO-IN)
      IVOT  = Q(JIN+2)
      JVOT  = LQ(JVOLUM-IVOT)
      IROTT = Q(JIN+4)
*
      IF (NLEVEL.GE.NLDEV(NLEVEL)) THEN
*       (case with JVOLUM structure locally developed)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         DO 169 ILEV = NLDEV(NLEVEL), NLEVEL
            IF (IQ(JPAR+1).EQ.0) THEN
               IF (ILEV.EQ.NLEVEL) THEN
                  JPAR = LQ(JPAR-IN)
               ELSE
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
               ENDIF
               IF (JPAR.EQ.0) GO TO 170
            ELSE IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
  169    CONTINUE
         JPAR = JPAR + 5
         GO TO 180
      ENDIF
*     (normal case)
  170 NPAR = Q(JVOT+5)
      IF (NPAR.EQ.0) THEN
         JPAR = JIN +9
      ELSE
         JPAR = JVOT +6
      ENDIF
*
*   * Compute distance to boundary of current content
*
  180 IF (IROTT.EQ.0) THEN
         DO 189 I = 1,3
            XT(I)   = XC(I) -Q(JIN+4+I)
            XT(I+3) = XC(I+3)
  189    CONTINUE
      ELSE
*       (later, code in line)
         CALL GITRAN (XC, Q(JIN+5), IROTT, XT)
         CALL GRMTD (XC(4), IROTT, XT(4))
      ENDIF
*
      IACT = 2
      ISHT = Q(JVOT+2)
      IF (ISHT.LT.5) THEN
         IF (ISHT.EQ.1) THEN
            CALL GNOBOX (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.2) THEN
            CALL GNOTRA(XT,Q(JPAR+1),IACT,1,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.3) THEN
            CALL GNOTRA(XT,Q(JPAR+1),IACT,2,SNEXT,SNXT,SAFE)
         ELSE
            CALL GNOTRP (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ENDIF
      ELSE IF (ISHT.LE.10) THEN
         IF (ISHT.EQ.5) THEN
            CALL GNOTUB(XT,Q(JPAR+1),IACT,1,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.6) THEN
            CALL GNOTUB(XT,Q(JPAR+1),IACT,2,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.7) THEN
            CALL GNOCON(XT,Q(JPAR+1),IACT,1,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.8) THEN
            CALL GNOCON(XT,Q(JPAR+1),IACT,2,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.9) THEN
            CALL GNOSPH (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ELSE
            CALL GNOPAR (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ENDIF
      ELSE IF (ISHT.EQ.11) THEN
         CALL GNOPGO (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.12) THEN
         CALL GNOPCO (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.13) THEN
         CALL GNOELT (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.14) THEN
         CALL GNOHYP (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.28) THEN
         CALL GSNGTR (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE,0)
      ELSE IF (ISHT.EQ.NSCTUB) THEN
         CALL GNOCTU (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE
         PRINT *, ' GNEXT : No code for shape ', ISHT
         STOP
      ENDIF
*
      IF (SAFE.LT.SAFETY) SAFETY = SAFE
      IF (SNXT.LT.SNEXT) THEN
         SNEXT = SNXT
         IF (ISEARC.EQ.-2) THEN
            IF (MOD(IQ(JSC0),2).NE.0) THEN
               IDIVB = IDIV
            ELSE
               DO 191 I = 1,3
                  X0(I) = XC(I) + SNXT*XC(I+3)
  191          CONTINUE
               IF (IAX.LE.3) THEN
                  IDIVB = LOCATF (Q(JSB+3), NSB, X0(IAX))
               ELSE
                  CALL GFCOOR (X0, IAX, CX)
                  IDIVB = LOCATF (Q(JSB+3), NSB, CX)
               ENDIF
               IF (IDIVB.LT.0) IDIVB = -IDIVB
               IF (IDIVB.EQ.0) THEN
                  IF (IAX.EQ.6) THEN
                     IDIVB = NSB
                  ELSE
                     IDIVB = 1
                  ENDIF
               ELSE IF (IDIVB.EQ.NSB) THEN
                  IF (IAX.NE.6) IDIVB = NSB - 1
               ENDIF
            ENDIF
         ENDIF
      ENDIF
*
      IF (ISEARC.EQ.-2) THEN
         IF (ICONT.EQ.NCONT) THEN
            IF (IDIVL.EQ.0) THEN
               IF (IDIVB.NE.0) THEN
                  IF (IDIV.EQ.IDIVB) GO TO 300
                  IF (.NOT.BTEST(IQ(JVO),2)) THEN
                     IDIVL = IDIVB
                     GO TO 193
                  ENDIF
               ENDIF
*
*   *         Compute distance to boundary of current volume
*
               JPAR = LQ(JGPAR-NLEVEL)
               IACT = 2
               ISH  = Q(JVO+2)
               IF (ISH.LT.5) THEN
                  IF (ISH.EQ.1) THEN
                     CALL GNBOX (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
                  ELSE IF (ISH.EQ.2) THEN
                     CALL GNTRAP (XC, Q(JPAR+1),IACT,1, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.3) THEN
                     CALL GNTRAP (XC, Q(JPAR+1),IACT,2, SNEXT,SNXT,SAFE)
                  ELSE
                     CALL GNTRP (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
                  ENDIF
               ELSE IF (ISH.LE.10) THEN
                  IF (ISH.EQ.5) THEN
                     CALL GNTUBE (XC, Q(JPAR+1),IACT,1, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.6) THEN
                     CALL GNTUBE (XC, Q(JPAR+1),IACT,2, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.7) THEN
                     CALL GNCONE (XC, Q(JPAR+1),IACT,1, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.8) THEN
                     CALL GNCONE (XC, Q(JPAR+1),IACT,2, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.9) THEN
                     CALL GNSPHR (XC, Q(JPAR+1),IACT, SNEXT, SNXT, SAFE)
                  ELSE
                     CALL GNPARA (XC, Q(JPAR+1),IACT, SNEXT, SNXT, SAFE)
                  ENDIF
               ELSE IF (ISH.EQ.12) THEN
                  CALL GNPCON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.11) THEN
                  CALL GNPGON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.13) THEN
                  CALL GNELTU (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.14) THEN
                  CALL GNHYPE (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.28) THEN
                  CALL GSNGTR (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE,1)
               ELSE IF (ISH.EQ.NSCTUB) THEN
                  CALL GNCTUB (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE
                  PRINT *, ' GNEXT : No code for shape ', ISH
                  STOP
               ENDIF
*
               IF (SAFE.LT.SAFETY) SAFETY = SAFE
               IF (SNXT.LT.SNEXT)  SNEXT  = SNXT
*
*   *         Check wether other pseudo-divisions have to be scanned
*
               DO 192 I = 1,3
                  X0(I) = XC(I) +SNXT*XC(I+3)
  192          CONTINUE
               IF (IAX.LE.3) THEN
                  IDIVL = LOCATF (Q(JSB+3), NSB, X0(IAX))
               ELSE
                  CALL GFCOOR (X0, IAX, CX)
                  IDIVL = LOCATF (Q(JSB+3), NSB, CX)
               ENDIF
               IF (IDIVL.LT.0) IDIVL = -IDIVL
               IF (IDIVL.EQ.0) THEN
                  IF(IAX.EQ.6)THEN
                     IDIVL=NSB
                  ELSE
                     IDIVL=1
                  ENDIF
               ELSEIF (IDIVL.EQ.NSB)THEN
                  IF(IAX.NE.6)IDIVL=NSB-1
               ENDIF
            ELSE
               IF (IDIV.EQ.IDIVB)   GO TO 400
            ENDIF
  193       IF ((IDIV-IDIVL)*INC.GE.0) GO TO 400
            IDIV = IDIV +INC
            GO TO 110
         ELSE
            ICONT = ICONT +1
            GO TO 140
         ENDIF
      ELSE
         IF (INEAR.EQ.NNEAR) GO TO 300
         INEAR = INEAR +1
         GO TO 130
      ENDIF
*
* ***    Case of volume incompletely divided
*
  200 JDIV  = LQ(JVO-1)
      IAXIS = Q(JDIV+1)
      IVOT  = Q(JDIV+2)
      JVOT  = LQ(JVOLUM-IVOT)
      ISHT  = Q(JVOT+2)
*
*  ** Get the division parameters
*
      IF (NLEVEL.LT.NLDEV(NLEVEL)) THEN
         JPARM = 0
      ELSE
*        (case with JVOLUM structure locally developed)
         JPARM = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         IF (NLEVEL.EQ.NLDEV(NLEVEL)) GO TO 215
         DO 210 ILEV = NLDEV(NLEVEL), NLEVEL-1
            IF (IQ(JPARM+1).EQ.0) THEN
               JPARM = LQ(JPARM-LINDEX(ILEV+1))
               IF (JPARM.EQ.0) GO TO 215
            ELSE IF (IQ(JPARM-3).GT.1) THEN
               JPARM = LQ(JPARM-LINDEX(ILEV+1))
            ELSE
               JPARM = LQ(JPARM-1)
            ENDIF
            IF (ILEV.EQ.NLEVEL-1) THEN
               NDIV = IQ(JPARM+1)
               ORIG =  Q(JPARM+2)
               SDIV =  Q(JPARM+3)
            ENDIF
  210    CONTINUE
         GO TO 220
      ENDIF
*     (normal case)
  215 NDIV = Q(JDIV+3)
      ORIG = Q(JDIV+4)
      SDIV = Q(JDIV+5)
*
*  ** Look at the first and the last divisions only
*
  220 IDT  = IDTYP(IAXIS, ISH)
      IF (IDT.EQ.1) THEN
         IN2 = 0
         IF (XC(IAXIS).LT.ORIG) THEN
            IN  = 1
         ELSE
            IN  = NDIV
         ENDIF
      ELSE IF (IDT.EQ.2) THEN
         R   = XC(1)**2 + XC(2)**2
         IF (ISH.EQ.9) R = R + XC(3)**2
         R   = SQRT(R)
         IN2 = 0
         IF (ISH.EQ.5.OR.ISH.EQ.6.OR.ISH.EQ.9) THEN
            IF (R.LT.ORIG) THEN
               IN  = 1
            ELSE
               IN  = NDIV
            ENDIF
         ELSE
            PRINT *, ' GNEXT : Partially divided ',ISH,IAXIS
            IN  = 1
            IF (NDIV.GT.1) IN2 = NDIV
         ENDIF
      ELSE IF (IDT.EQ.4) THEN
         IN2 = 0
         RXY = XC(1)**2 + XC(2)**2
         RXY = SQRT(RXY)
         IF (XC(3).NE.0.0) THEN
            THET = RADDEG * ATAN (RXY/XC(3))
            IF (THET.LT.0.0) THET = THET + 180.0
         ELSE
            THET = 90.
         ENDIF
         IF (THET.LE.ORIG) THEN
            IN  = 1
         ELSE
            IN  = NDIV
         ENDIF
      ELSE
         PRINT *, ' GNEXT : Partially divided ',ISH,IAXIS
         IN2 = 0
         IF (ISH.EQ.5.OR.ISH.EQ.7) THEN
            IN  = 1
            IF (NDIV.GT.1) IN2 = NDIV
         ELSE
            IF (XC(1).NE.0.0.OR.XC(2).NE.0.0) THEN
               PHI = RADDEG * ATAN2 (XC(2), XC(1))
            ELSE
               PHI = 0.0
            ENDIF
            IF (ISH.EQ.6.OR.ISH.EQ.8) THEN
               IF (PHI.LT.ORIG) THEN
                  IN  = 1
               ELSE
                  IN  = NDIV
               ENDIF
            ELSE
               IN  = 1
               IF (NDIV.GT.1) IN2 = NDIV
            ENDIF
         ENDIF
      ENDIF
*
  225 IF (IDT.EQ.1) THEN
         DO 231 I = 1, 3
            X0(I) = 0.0
  231    CONTINUE
         X0(IAXIS) = ORIG + (IN - 0.5) * SDIV
         IF (ISH.EQ.4.OR.(ISH.EQ.10.AND.IAXIS.NE.1)) THEN
            CALL GCENT (IAXIS, X0)
         ENDIF
         DO 232 I = 1, 3
            XT(I)   = XC(I) - X0(I)
            XT(I+3) = XC(I+3)
  232    CONTINUE
      ELSE IF (IDT.EQ.3) THEN
         PH0  = DEGRAD * (ORIG + (IN - 0.5) * SDIV)
         CPHR = COS(PH0)
         SPHR = SIN(PH0)
         DO 233 I = 1, 4, 3
            XT(I)   = XC(I)*CPHR + XC(I+1)*SPHR
            XT(I+1) = XC(I+1)*CPHR - XC(I)*SPHR
            XT(I+2) = XC(I+2)
  233    CONTINUE
      ELSE
         DO 234 I = 1, 6
            XT(I) = XC(I)
  234    CONTINUE
      ENDIF
*
      IF (JPARM.NE.0) THEN
         IF (IQ(JPARM-3).GT.1) THEN
            JPAR = LQ(JPARM-IN)
         ELSE
            JPAR = LQ(JPARM-1)
         ENDIF
         JPAR = JPAR + 5
      ELSE
         JPAR = JVOT + 6
      ENDIF
*
      IACT = 2
      IF (ISHT.LT.5) THEN
         IF (ISHT.EQ.1) THEN
            CALL GNOBOX (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.2) THEN
            CALL GNOTRA (XT, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.3) THEN
            CALL GNOTRA (XT, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNOTRP (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISHT.LE.10) THEN
         IF (ISHT.EQ.5) THEN
            CALL GNOTUB (XT, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.6) THEN
            CALL GNOTUB (XT, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.7) THEN
            CALL GNOCON (XT, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.8) THEN
            CALL GNOCON (XT, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.9) THEN
            CALL GNOSPH (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNOPAR (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISHT.EQ.11) THEN
         CALL GNOPGO (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.12) THEN
         CALL GNOPCO (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.13) THEN
         CALL GNOELT (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.14) THEN
         CALL GNOHYP (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.28) THEN
         CALL GSNGTR (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE,0)
      ELSE IF (ISHT.EQ.NSCTUB) THEN
         CALL GNOCTU (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE
         PRINT *, ' GNEXT : No code for shape ', ISHT
         STOP
      ENDIF
*
      IF (SAFE.LT.SAFETY) SAFETY = SAFE
      IF (SNXT.LT.SNEXT)  SNEXT  = SNXT
*
      IF (IN2.NE.0) THEN
         IF (IN2.NE.IN) THEN
            IN  = IN2
            GO TO 225
         ENDIF
      ENDIF
*
* ***  Calculate SNEXT and SAFETY with respect to the Mother
* ***            SAFETY only for concave volumes if finite SNEXT
* ***            has been found with respect to one of its contents
*
  300 IACT = 2
      IF (SNEXT.LT.0.9*BIG) THEN
         IF (.NOT.BTEST(IQ(JVO),2)) IACT = 0
      ENDIF
      JPAR = LQ(JGPAR-NLEVEL)
      IF (ISH.LT.5) THEN
         IF (ISH.EQ.1) THEN
            CALL GNBOX (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE )
         ELSE IF (ISH.EQ.2) THEN
            CALL GNTRAP (XC, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.3) THEN
            CALL GNTRAP (XC, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNTRP (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISH.LE.10) THEN
         IF (ISH.EQ.5) THEN
            CALL GNTUBE (XC, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.6) THEN
            CALL GNTUBE (XC, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.7) THEN
            CALL GNCONE (XC, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.8) THEN
            CALL GNCONE (XC, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.9) THEN
            CALL GNSPHR (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNPARA (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISH.EQ.12) THEN
         CALL GNPCON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.11) THEN
         CALL GNPGON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.13) THEN
         CALL GNELTU (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.14) THEN
         CALL GNHYPE (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.28) THEN
         CALL GSNGTR (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE,1)
      ELSE IF (ISH.EQ.NSCTUB) THEN
         CALL GNCTUB (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE
         PRINT *, ' GNEXT : No code for shape ', ISH
         STOP
      ENDIF
*
      IF (SAFE.LT.SAFETY) SAFETY = SAFE
      IF (SNXT.LT.SNEXT)  SNEXT  = SNXT
*
  400 IF (GONLY(NLEVEL).EQ.0.) THEN
*
* ***   Case of a 'NOT ONLY' volume -> step search
*
         SAFETY = 0.
         EPSI2  = 0.5*EPSIL
         ST     = SNEXT -EPSI2
         IF (ST.LE.0) GO TO 999
         EPSI3  = 10.*EPSIL
         IF (ST.LE.EPSI3) THEN
            NN = 1
         ELSE
            NN = ST/EPSI3 +1
            ST = ST/NN
         ENDIF
*
         NBIN = 0
         SN   = 0.
  420    SN   = SN +ST
         DO 429 I = 1,3
            XT(I) = X(I) +SN*X(I+3)
  429    CONTINUE
*
         CALL GINVOL (XT, ISAME)
         IF (ISAME.EQ.0) THEN
            IF (ST.LT.EPSI2) GO TO 490
            SN   = SN -ST
            ST   = 0.5*ST
            NBIN = 1
            GO TO 420
         ENDIF
*
         IF (NBIN.NE.0) THEN
            IF (ST.LT.EPSI2) THEN
               SN = SN +EPSI2
               GO TO 490
            ELSE
               ST = 0.5*ST
               GO TO 420
            ENDIF
         ENDIF
         NN = NN -1
         IF (NN.GT.0) GO TO 420
         GO TO 999
*
  490    IF (SN.LT.SNEXT) SNEXT = SN
      ENDIF
*                                                              END GNEXT
  999 END
 
#endif

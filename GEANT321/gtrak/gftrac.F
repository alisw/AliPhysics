*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GFTRAC
C.    ******************************************************************
C.    *                                                                *
C.    *   SUBR. GFTRAC                                                 *
C.    *                                                                *
C.    *   Selects next track segment to be processed and extracts from *
C.    *    the stack JTRACK the relevant information to reload commons *
C.    *                                                                *
C.    *   Called by : GTREVE                                           *
C.    *   Authors   : S.Banerjee, F.Bruyant                            *
C.    *                                                                *
C.    ******************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gcnum.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcpoly.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      REAL      XC(3), XT(3), X0(3)
      INTEGER   IDTYP(3,12)
      LOGICAL   BTEST
C.
      SAVE MANY
      DATA  MANY / 0/
      DATA  IDTYP / 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 2, 3, 1,
     +              2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 4, 3, 1, 1, 1,
     +              2, 3, 1, 2, 3, 1/
C.    ------------------------------------------------------------------
*
* *** Process next track in 'IN current VOlume' chain, if any
*
      IF (NJTMAX.LT.0) THEN
*
*  **   Reactivate parallel tracking if enough space available
*
         IF (NALIVE.LE.NJTMIN) NJTMAX = -NJTMAX
*
*  **   Update common /GCVOLU/ and structure JGPAR if necessary
*
         NLEVEL = NLVSAV
         ISKP   = 1
         DO 9 ILEV = 2,NLDOWN
            IF (ISKP.NE.0) THEN
               IF (LINDEX(ILEV).EQ.LINSAV(ILEV)) GO TO 9
               ISKP = 0
            ENDIF
            JSKLD = LQ(JSKLT-ILEV)
            JSKD  = LQ(JSKLD-LINSAV(ILEV))
            IVO   = IQ(JSKD+2)
            LQ(JGPAR-ILEV) = LQ(JSKD-1)
            IQ(JGPAR+ILEV) = IQ(JSKD+1)
            LVOLUM(ILEV)   = IVO
            NAMES(ILEV)    = IQ(JVOLUM+IVO)
            LINDEX(ILEV)   = LINSAV(ILEV)
            LINMX(ILEV)    = LMXSAV(ILEV)
            JVOM = LQ(JVOLUM-LVOLUM(ILEV-1))
            IF (Q(JVOM+3).GT.0.) THEN
               JIN = LQ(JVOM-LINDEX(ILEV))
               NUMBER(ILEV) = Q(JIN+3)
               GONLY(ILEV)  = Q(JIN+8)
            ELSE
               NUMBER(ILEV) = LINDEX(ILEV)
               GONLY(ILEV)  = GONLY(ILEV-1)
            ENDIF
            IF (LQ(LQ(JVOLUM-IVO)).EQ.0) THEN
               NLDEV(ILEV) = NLDEV(ILEV-1)
            ELSE
               NLDEV(ILEV) = ILEV
            ENDIF
            GTRAN(1,ILEV) = Q(JSKD+3)
            GTRAN(2,ILEV) = Q(JSKD+4)
            GTRAN(3,ILEV) = Q(JSKD+5)
            DO 8 I = 1, 10, 2
               GRMAT(I,ILEV) = Q(JSKD+5+I)
               GRMAT(I+1,ILEV) = Q(JSKD+6+I)
    8       CONTINUE
    9    CONTINUE
* **
         IF (NJINVO.NE.0) GO TO 800
         IFUPD = 0
      ELSE
         IF (NJINVO.NE.0) GO TO 800
         IFUPD = 1
      ENDIF
*
* *** 'IN current VOlume' chain is empty, refill from JSKLT structure
*     Scan brother chains, starting from current one when going up in
*     the skeleton structure
*
   10 INSK  = 1
*
   11 NLEVEL = NLDOWN
      JSKLD  = LQ(JSKLT-NLEVEL)
      NINSK  = LINMX(NLEVEL)
      IDO    = 1
*
   20 IF (IQ(JSKLD+INSK).EQ.0) GO TO 589
      JSKD = LQ(JSKLD-INSK)
      IVO  = IQ(JSKD+2)
      IF (IFUPD.NE.0.AND.NLEVEL.GT.1) THEN
*
*  **   Update common /GCVOLU/ for level NLEVEL
*
         LQ(JGPAR-NLEVEL) = LQ(JSKD-1)
         IQ(JGPAR+NLEVEL) = IQ(JSKD+1)
         LVOLUM(NLEVEL)   = IVO
         NAMES(NLEVEL)    = IQ(JVOLUM+IVO)
         LINDEX(NLEVEL)   = INSK
         JVOM = LQ(JVOLUM-LVOLUM(NLEVEL-1))
         IF (Q(JVOM+3).GT.0.) THEN
            JIN = LQ(JVOM-INSK)
            NUMBER(NLEVEL) = Q(JIN+3)
            GONLY(NLEVEL)  = Q(JIN+8)
         ELSE
            NUMBER(NLEVEL) = INSK
            GONLY(NLEVEL)  = GONLY(NLEVEL-1)
         ENDIF
         IF (LQ(LQ(JVOLUM-IVO)).EQ.0) THEN
            NLDEV(NLEVEL) = NLDEV(NLEVEL-1)
         ELSE
            NLDEV(NLEVEL) = NLEVEL
         ENDIF
         GTRAN(1,NLEVEL) = Q(JSKD+3)
         GTRAN(2,NLEVEL) = Q(JSKD+4)
         GTRAN(3,NLEVEL) = Q(JSKD+5)
         DO 29 I = 1, 10, 2
            GRMAT(I,NLEVEL) = Q(JSKD+5+I)
            GRMAT(I+1,NLEVEL) = Q(JSKD+6+I)
   29    CONTINUE
      ENDIF
*
      JVO  = LQ(JVOLUM-IVO)
      IF (Q(JVO+3).EQ.0.) GO TO 600
      NIN = Q(JVO+3)
*
*  ** Sort-out unsorted-out elements in first non-empty brother chain
*
      LPREV = JSKLD +INSK
      NCUR  = IQ(LPREV)
   50 LCUR = JTRACK +(NCUR-1)*NWTRAC
      IF (IQ(LCUR+2).NE.0) GO TO 600
      NSTO = IQ(LCUR+1)
*
      IPCUR = LCUR +NWINT
C*****  Code Expanded From Routine:  GTRNSF
C
      IF (GRMAT(10,NLEVEL) .EQ. 0.) THEN
         XC(1) = Q(1+IPCUR) - GTRAN(1,NLEVEL)
         XC(2) = Q(2+IPCUR) - GTRAN(2,NLEVEL)
         XC(3) = Q(3+IPCUR) - GTRAN(3,NLEVEL)
*
      ELSE
         XL11X = Q(1+IPCUR) - GTRAN(1,NLEVEL)
         XL21X = Q(2+IPCUR) - GTRAN(2,NLEVEL)
         XL31X = Q(3+IPCUR) - GTRAN(3,NLEVEL)
         XC(1) = XL11X*GRMAT(1,NLEVEL) + XL21X*GRMAT(2,NLEVEL) + XL31X*
     1      GRMAT(3,NLEVEL)
         XC(2) = XL11X*GRMAT(4,NLEVEL) + XL21X*GRMAT(5,NLEVEL) + XL31X*
     1      GRMAT(6,NLEVEL)
         XC(3) = XL11X*GRMAT(7,NLEVEL) + XL21X*GRMAT(8,NLEVEL) + XL31X*
     1      GRMAT(9,NLEVEL)
 
      ENDIF
C*****  End of Code Expanded From Routine:  GTRNSF
*
      IF (NIN.LT.0) GO TO 200
*
*   *   Case with contents defined by Position
*
      JNEAR = LQ(JVO-NIN-1)
      INFROM = IQ(LCUR+11)
      IF (INFROM.GT.0) THEN
         JIN   = LQ(JVO-INFROM)
         IF (LQ(JIN-1).NE.0) JNEAR = LQ(JIN-1)
      ENDIF
      IF (IQ(JNEAR+2).EQ.0) GO TO 300
      ISEARC = Q(JVO+1)
      IF (ISEARC.LT.0) THEN
*
*       Prepare access list when contents have been ordered by GSORD
*
         JSB = LQ(LQ(JVO-NIN-1))
         IAX = Q(JSB+1)
         NSB = Q(JSB+2)
         IF (IAX.LE.3) THEN
            IDIV = LOCATF (Q(JSB+3), NSB, XC(IAX))
         ELSE
            CALL GFCOOR (XC, IAX, CX)
            IDIV = LOCATF (Q(JSB+3), NSB, CX)
         ENDIF
         IF (IDIV.LT.0) IDIV = -IDIV
         IF (IDIV.EQ.0) THEN
            IF (IAX.NE.6) GO TO 300
            IDIV = NSB
         ELSE IF (IDIV.EQ.NSB) THEN
            IF (IAX.NE.6) GO TO 300
         ENDIF
         JSC0 = LQ(JVO-NIN-2)
         NCONT = IQ(JSC0+IDIV)
         IF (NCONT.LE.0) GO TO 300
         JSCV = LQ(JSC0-IDIV)
         ICONT = 1
         GO TO 120
      ELSE
         IF (ISEARC.GT.0) THEN
#if !defined(CERNLIB_USRJMP)
            CALL GUNEAR (ISEARC, 1, XC, JNEAR)
#endif
#if defined(CERNLIB_USRJMP)
            CALL JUMPT4(JUNEAR, ISEARC, 1, XC, JNEAR)
#endif
            IF (IQ(JNEAR+1).EQ.0) GO TO 300
         ENDIF
         JNEAR  = JNEAR +1
         NNEAR  = IQ(JNEAR)
         INEAR  = 1
      ENDIF
*
  110 IN = IQ(JNEAR+INEAR)
      IF (IN.GT.0) GO TO 150
      GO TO 190
*
  120 IN = IQ(JSCV+ICONT)
*
*     For each selected content in turn, check if point is in
*
  150 JIN  = LQ(JVO-IN)
      IVOT = Q(JIN+2)
      JVOT = LQ(JVOLUM-IVOT)
      IF (BTEST(IQ(JVOT),1)) THEN
*       (case with JVOLUM structure locally developed)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         DO 169 ILEV = NLDEV(NLEVEL), NLEVEL
            IF (IQ(JPAR+1).EQ.0) THEN
               IF (ILEV.EQ.NLEVEL) THEN
                  JPAR = LQ(JPAR-IN)
               ELSE
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
               ENDIF
            ELSE IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
  169    CONTINUE
         JPAR = JPAR +5
         NPAR = IQ(JPAR)
         GO TO 175
      ENDIF
*      (normal case)
      NPAR = Q(JVOT+5)
      IF (NPAR.EQ.0) THEN
         JPAR = JIN +9
         NPAR = Q(JPAR)
      ELSE
         JPAR = JVOT +6
      ENDIF
*
  175 IROTT  = Q(JIN+4)
C*****  Code Expanded From Routine:  GITRAN
C.
C.    ------------------------------------------------------------------
C.
      IF (IROTT.EQ.0) THEN
         XT(1) = XC(1) - Q(JIN+5)
         XT(2) = XC(2) - Q(JIN+6)
         XT(3) = XC(3) - Q(JIN+7)
*
      ELSE
         XL1 = XC(1) - Q(5+JIN)
         XL2 = XC(2) - Q(6+JIN)
         XL3 = XC(3) - Q(7+JIN)
         JR = LQ(JROTM-IROTT)
         XT(1) = XL1*Q(JR+1) + XL2*Q(JR+2) + XL3*Q(JR+3)
         XT(2) = XL1*Q(JR+4) + XL2*Q(JR+5) + XL3*Q(JR+6)
         XT(3) = XL1*Q(JR+7) + XL2*Q(JR+8) + XL3*Q(JR+9)
*
      ENDIF
C*****  Code Expanded From Routine:  GITRAN
      CALL GINME (XT, Q(JVOT+2), Q(JPAR+1), IYES)
      IF (IYES.NE.0) THEN
*
*       Volume found at deeper level
*
         NLDOWN = NLEVEL +1
         LINMX(NLDOWN) = NIN
         JSKL = LQ(JSKLT-NLDOWN)
*
*       Clear skeleton at lowest level if necessary
*
         JOFF = JSKL +IQ(JSKL-3)
         DO 184 ILEV = 1,NLEVEL
            IF (IQ(JOFF+ILEV).EQ.LINDEX(ILEV)) GO TO 184
            DO 182 I = ILEV,NLEVEL
               IQ(JOFF+I) = LINDEX(I)
  182       CONTINUE
            DO 183 I = 1,NIN
               JSK = LQ(JSKL-I)
               IQ(JSK+1) = 0
  183       CONTINUE
            GO TO 185
  184    CONTINUE
*
*       Prepare skeleton for level down if not yet done
*
  185    JSK  = LQ(JSKL-IN)
         IF (IQ(JSK+1).EQ.0) THEN
            LQ(JSK-1) = JPAR
            IQ(JSK+1) = NPAR
            IQ(JSK+2) = IVOT
            CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL),
     +                   Q(JIN+5), IROTT, Q(JSK+3), Q(JSK+6))
         ENDIF
         GO TO 500
*
      ENDIF
*
  190 IF (ISEARC.LT.0) THEN
         IF (ICONT.EQ.NCONT) GO TO 300
         ICONT = ICONT +1
         GO TO 120
      ELSE
         IF (INEAR.EQ.NNEAR) GO TO 300
         INEAR = INEAR +1
         GO TO 110
      ENDIF
*
*   *   Case with contents defined by division
*
  200 JDIV  = LQ(JVO-1)
      ISH   = Q(JVO+2)
      IAXIS = Q(JDIV+1)
      IVOT  = Q(JDIV+2)
      JVOT  = LQ(JVOLUM-IVOT)
      IF (NLEVEL.LT.NLDEV(NLEVEL)) THEN
         JPAR = 0
      ELSE
*       (case with structure JVOLUM locally developped)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         IF (NLEVEL.EQ.NLDEV(NLEVEL)) GO TO 250
         DO 249 ILEV = NLDEV(NLEVEL), NLEVEL-1
            IF (IQ(JPAR+1).EQ.0) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
               IF (JPAR.EQ.0) GO TO 250
            ELSE IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
            IF (ILEV.EQ.NLEVEL-1) THEN
               NDIV  = IQ(JPAR+1)
               ORIG  =  Q(JPAR+2)
               SDIV  =  Q(JPAR+3)
            ENDIF
  249    CONTINUE
         GO TO 260
      ENDIF
*      (normal case)
  250 NDIV  = Q(JDIV+3)
      ORIG  = Q(JDIV+4)
      SDIV  = Q(JDIV+5)
*
  260 IDT = IDTYP(IAXIS,ISH)
      IF (IDT.EQ.1) THEN
*
*      Division along X, Y or Z axis
*
         XTT = XC(IAXIS)
         IF (ISH.EQ.10) THEN
            IF (IAXIS.NE.3) THEN
               XTT = XTT - Q(LQ(JGPAR-NLEVEL)+IAXIS+4) * XC(3)
               IF (IAXIS.EQ.1) THEN
                  YT  = XC(2) - Q(LQ(JGPAR-NLEVEL)+6) * XC(3)
                  XTT = XTT - Q(LQ(JGPAR-NLEVEL)+4) * YT
               ENDIF
            ENDIF
         ENDIF
         IN = (XTT -ORIG)/SDIV +1
      ELSE IF (IDT.EQ.2) THEN
*
*       Division along R axis
*
         R = XC(1)**2 + XC(2)**2
         IF (ISH.EQ.9) R = R + XC(3)**2
         R = SQRT (R)
         IF (ISH.EQ.5.OR.ISH.EQ.6.OR.ISH.EQ.9) THEN
            IN = (R - ORIG) / SDIV + 1
         ELSE IF (ISH.EQ.7.OR.ISH.EQ.8) THEN
            IPAR = LQ(JGPAR-NLEVEL)
            DR   = 0.5 * (Q(IPAR+4) - Q(IPAR+2)) / Q(IPAR+1)
            RMN  = 0.5 * (Q(IPAR+4) + Q(IPAR+2)) + DR * XC(3)
            DR   = 0.5 * (Q(IPAR+5) - Q(IPAR+3)) / Q(IPAR+1)
            RMX  = 0.5 * (Q(IPAR+5) + Q(IPAR+3)) + DR * XC(3)
            STP  = (RMX - RMN) / NDIV
            IN   = (R - RMN) / STP + 1
         ELSE
            IPAR = LQ(JGPAR-NLEVEL)
            IF (ISH.EQ.12) THEN
               IPT = IPAR + 1
            ELSE
               IPT = IPAR + 2
            ENDIF
            IF (IZSEC.GT.0) THEN
               IPT = IPT + 3 * IZSEC
            ELSE
               NZ  = Q(IPT+2)
               DO 261 IZ = 1, NZ-1
                  IF((XC(3)-Q(IPT+3*IZ))*(XC(3)-Q(IPT+3*IZ+3)).LE.0.)
     +            THEN
                     IZSEC = IZ
                     IPT   = IPT + 3 * IZSEC
                     GO TO 262
                  ENDIF
  261          CONTINUE
               IN  = 0
               GO TO 265
            ENDIF
  262       POR1 = (Q(IPT+3) - XC(3)) / (Q(IPT+3) - Q(IPT))
            POR2 = (XC(3) - Q(IPT)) / (Q(IPT+3) - Q(IPT))
            RMN  = Q(IPT+1) * POR1 + Q(IPT+4) * POR2
            RMX  = Q(IPT+2) * POR1 + Q(IPT+5) * POR2
            IF (ISH.EQ.11) THEN
               NPDV = Q(IPAR+3)
               DPH  = Q(IPAR+2) / NPDV
               IF (IPSEC.LE.0) THEN
                  IF (XC(1).NE.0..OR.XC(2).NE.0.) THEN
                     PHI  = RADDEG * ATAN2 (XC(2), XC(1))
                  ELSE
                     PHI  = 0.0
                  ENDIF
                  PH0  = MOD (PHI-Q(IPAR+1)+360., 360.)
                  IPSEC= PH0/DPH + 1
               ENDIF
               PH   = DEGRAD * (Q(IPAR+1) + (IPSEC - 0.5) * DPH)
               R    = XC(1) * COS(PH) + XC(2) * SIN(PH)
            ENDIF
            STP = (RMX - RMN) / NDIV
            IN  = (R - RMN) / STP + 1
         ENDIF
      ELSE IF (IDT.EQ.3) THEN
*
*       Division along Phi axis
*
         IF (XC(1).NE.0..OR.XC(2).NE.0.) THEN
            PHI = RADDEG * ATAN2 (XC(2), XC(1))
         ELSE
            PHI = 0.
         ENDIF
         IN  = MOD (PHI-ORIG+360., 360.) / SDIV + 1
      ELSE IF (IDT.EQ.4) THEN
*
*       Division along Theta axis
*
         IF (XC(3).NE.0.0) THEN
            RXY  = SQRT (XC(1)**2 + XC(2)**2)
            THET = RADDEG * ATAN (RXY/XC(3))
            IF (THET.LT.0.0)  THET = THET + 180.0
         ELSE
            THET = 90.0
         ENDIF
         IN   = (THET - ORIG) / SDIV + 1
      ENDIF
*
  265 IF (IN.GT.NDIV) IN = 0
      IF (IN.LE.0) GO TO 300
*
      IF (JPAR.NE.0) THEN
         IF (IQ(JPAR-3).GT.1) THEN
            JPAR = LQ(JPAR-IN)
         ELSE
            JPAR = LQ(JPAR-1)
         ENDIF
         JPAR = JPAR + 5
         NPAR = IQ(JPAR)
      ELSE
         NPAR = Q(JVOT+5)
         JPAR = JVOT + 6
      ENDIF
*
*      Volume found at deeper level
*
      NLDOWN = NLEVEL +1
      LINMX(NLDOWN) = NDIV
      JSKL = LQ(JSKLT-NLDOWN)
*
*      Clear skeleton at lowest level if necessary
*
      JOFF = JSKL +IQ(JSKL-3)
      DO 269 ILEV = 1,NLEVEL
         IF (IQ(JOFF+ILEV).EQ.LINDEX(ILEV)) GO TO 269
         DO 267 I = ILEV,NLEVEL
            IQ(JOFF+I) = LINDEX(I)
  267    CONTINUE
         DO 268 I = 1,NDIV
            JSK = LQ(JSKL-I)
            IQ(JSK+1) = 0
  268    CONTINUE
         GO TO 270
  269 CONTINUE
*
*       Prepare skeleton at level down if not yet done
*
  270 JSK  = LQ(JSKL-IN)
      IF (IQ(JSK+1).EQ.0) THEN
         LQ(JSK-1) = JPAR
         IQ(JSK+1) = NPAR
         IQ(JSK+2) = IVOT
*
         IF (IDT.EQ.1) THEN
            X0(1) = 0.0
            X0(2) = 0.0
            X0(3) = 0.0
            X0(IAXIS) = ORIG + (IN - 0.5) * SDIV
            IF (ISH.EQ.4.OR.(ISH.EQ.10.AND.IAXIS.NE.1)) THEN
               CALL GCENT (IAXIS, X0)
            ENDIF
            IF (GRMAT(10,NLEVEL).EQ.0.0) THEN
               Q(JSK+3) = GTRAN(1,NLEVEL) + X0(1)
               Q(JSK+4) = GTRAN(2,NLEVEL) + X0(2)
               Q(JSK+5) = GTRAN(3,NLEVEL) + X0(3)
               DO 278 I = 1, 10, 2
                  Q(JSK+5+I) = GRMAT(I,NLEVEL)
                  Q(JSK+6+I) = GRMAT(I+1,NLEVEL)
  278          CONTINUE
            ELSE
               CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), X0, 0,
     +                      Q(JSK+3), Q(JSK+6))
            ENDIF
*
         ELSE IF (IDT.EQ.3.OR.IDT.EQ.4) THEN
            IF (IDT.EQ.3) THEN
               PH0  = DEGRAD * (ORIG + (IN - 0.5) * SDIV)
               CPHR = COS (PH0)
               SPHR = SIN (PH0)
            ELSE
               PH0  = 0.0
               CPHR = 1.0
               SPHR = 0.0
            ENDIF
            DO 279 I = 1, 3
               Q(JSK+2+I) = GTRAN(I,NLEVEL)
               Q(JSK+5+I) = GRMAT(I,NLEVEL)*CPHR +GRMAT(I+3,NLEVEL)*SPHR
               Q(JSK+8+I) = GRMAT(I+3,NLEVEL)*CPHR -GRMAT(I,NLEVEL)*SPHR
               Q(JSK+11+I)= GRMAT(I+6,NLEVEL)
  279       CONTINUE
            IF (PH0.EQ.0.0.AND.GRMAT(10,NLEVEL).EQ.0.0) THEN
               Q(JSK+15) = 0.0
            ELSE
               Q(JSK+15) = 1.0
            ENDIF
            IF (ISH.EQ.11) IPSEC = 1
*
         ELSE
            Q(JSK+3) = GTRAN(1,NLEVEL)
            Q(JSK+4) = GTRAN(2,NLEVEL)
            Q(JSK+5) = GTRAN(3,NLEVEL)
            DO 281 I = 1, 10, 2
               Q(JSK+5+I) = GRMAT(I,NLEVEL)
               Q(JSK+6+I) = GRMAT(I+1,NLEVEL)
  281       CONTINUE
         ENDIF
*
      ENDIF
      GO TO 500
*
  300 IF (GONLY(NLEVEL).EQ.0.) THEN
         IF (MANY.EQ.0) THEN
            WRITE (CHMAIL, 1001)
            CALL GMAIL (0 ,0)
            MANY = 1
         ENDIF
      ENDIF
*
      IQ(LCUR+2) = 1
      LPREV = LCUR +1
      GO TO 510
*
*       Move track down in skeleton
*
  500 IQ(LPREV)   = NSTO
      IQ(LCUR+1)  = IQ(JSKL+IN)
*      (reset INFROM to 0)
      IQ(LCUR+11) = 0
      IQ(JSKL+IN) = NCUR
*
  510 IF (NSTO.EQ.0) THEN
         GO TO 600
      ELSE
         NCUR = NSTO
         GO TO 50
      ENDIF
*
  589 IF (IDO.LT.NINSK) THEN
         IDO  = IDO +1
         INSK = INSK +1
         IF (INSK.GT.NINSK) INSK = 1
         IFUPD = 1
         GO TO 20
      ENDIF
*
*  **  No more elements at lowest level, go one level up in skeleton
*
      NLDOWN = NLDOWN -1
      INSK   = LINDEX(NLDOWN)
      IFUPD  = 0
      GO TO 11
*
  600 IF (NLDOWN.GT.NLEVEL) THEN
         IFUPD = 1
         GO TO 10
      ENDIF
*
*  **  Prepare 'IN current VOlume' chain
*
      NJINVO = IQ(JSKLD+INSK)
      IQ(JSKLD+INSK) = 0
*
      IF (NJTMAX.LT.0) THEN
*       (save status of skeleton for later reactivation of // tracking)
         DO 609 I = 2,NLEVEL
            LINSAV(I) = LINDEX(I)
            LMXSAV(I) = LINMX(I)
  609    CONTINUE
      ENDIF
*
* *** Fetch information for next track segment to be processed
*
  800 NCUR   = NJINVO
      LCUR   = JTRACK +(NCUR-1)*NWTRAC
      NJINVO = IQ(LCUR+1)
      NTMULT = IQ(LCUR+3)
      ITRA   = IQ(LCUR+4)
      ISTAK  = IQ(LCUR+5)
      IPART  = IQ(LCUR+6)
      NSTEP  = IQ(LCUR+7)
*free IDECAD = IQ(LCUR+8)
      IEKBIN = IQ(LCUR+9)
      ISTORY = IQ(LCUR+10)
      INFROM = IQ(LCUR+11)
*
      IF (IPART.NE.IPAOLD) THEN
         JPA = LQ(JPART-IPART)
         DO 819 I = 1,5
            NAPART(I) = IQ(JPA+I)
  819    CONTINUE
         ITRTYP = Q(JPA+6)
         AMASS  = Q(JPA+7)
         CHARGE = Q(JPA+8)
         TLIFE  = Q(JPA+9)
         IPAOLD = IPART
         IUPD   = 0
      ENDIF
*
      IPCUR = LCUR +NWINT
      DO 829 I = 1,7
         VECT(I) = Q(IPCUR+I)
  829 CONTINUE
      GEKIN  = Q(IPCUR+8)
      SLENG  = Q(IPCUR+9)
      GEKRAT = Q(IPCUR+10)
      TOFG   = Q(IPCUR+11)
      UPWGHT = Q(IPCUR+12)
*
      GETOT  = GEKIN +AMASS
      SAFETY = 0.
*
      IPCUR = IPCUR +NWREAL
      IF (ITRTYP.EQ.1) THEN
*      Photons
         ZINTPA = Q(IPCUR+1)
         ZINTCO = Q(IPCUR+2)
         ZINTPH = Q(IPCUR+3)
         ZINTPF = Q(IPCUR+4)
         ZINTRA = Q(IPCUR+5)
      ELSE IF (ITRTYP.EQ.2) THEN
*     Electrons
         ZINTBR = Q(IPCUR+1)
         ZINTDR = Q(IPCUR+2)
         ZINTAN = Q(IPCUR+3)
      ELSE IF (ITRTYP.EQ.3) THEN
*     Neutral hadrons
         SUMLIF = Q(IPCUR+1)
         ZINTHA = Q(IPCUR+2)
      ELSE IF (ITRTYP.EQ.4) THEN
*     Charged hadrons
         SUMLIF = Q(IPCUR+1)
         ZINTHA = Q(IPCUR+2)
         ZINTDR = Q(IPCUR+3)
      ELSE IF (ITRTYP.EQ.5) THEN
*     Muons
         SUMLIF = Q(IPCUR+1)
         ZINTBR = Q(IPCUR+2)
         ZINTPA = Q(IPCUR+3)
         ZINTDR = Q(IPCUR+4)
         ZINTMU = Q(IPCUR+5)
      ELSE IF (ITRTYP.EQ.7) THEN
*     Cerenkov photons
         ZINTLA = Q(IPCUR+1)
      ELSE IF (ITRTYP.EQ.8) THEN
*     Ions
         ZINTHA = Q(IPCUR+1)
         ZINTDR = Q(IPCUR+2)
      ENDIF
*
*   *  Reset NUMED
*
      JVO   = LQ(JVOLUM-LVOLUM(NLEVEL))
      NUMED = Q(JVO+4)
*
*     Link selected track segment area to 'garbaged' chain
*
      IQ(LCUR+1) = NJGARB
      NJGARB     = NCUR
*
*     Save skeleton status when parallel tracking is frozen
*
      IF (NJTMAX.LT.0) THEN
         NLVSAV = NLEVEL
         DO 889 ILEV = 2,NLDOWN
            LINSAV(ILEV) = LINDEX(ILEV)
            LMXSAV(ILEV) = LINMX(ILEV)
  889    CONTINUE
      ENDIF
*
 1001 FORMAT (' GFTRAC : Simple NOT-ONLY configuration assumed. OK?')
*                                                             END GFTRAC
      END

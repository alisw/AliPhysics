*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GMEDIV (JVO, IN, XC, IFL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Updates the common /GCVOLU/ and the structure JGPAR          *
C.    *     for contents defined by division.                          *
C.    *                                                                *
C.    *   For IFL nonzero, it also checks if the point XC is inside    *
C.    *     the content. It returns IN = 0, if the point is outside.   *
C.    *     Otherwise, it transforms XC in the local system.           *
C.    *                                                                *
C.    *   For IFL zero, IN is returned 0, if IN > NDIV.                *
C.    *                                                                *
C.    *   Input : JVO, IN, XC, IFL                                     *
C.    *   Output : IN, XC                                              *
C.    *                                                                *
C.    *   Called by : GDRAW, GMEDIA                                    *
C.    *   Authors   : S.Banerjee, R.Brun, F.Bruyant, A.McPherson       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcpoly.inc"
#include "geant321/gcvolu.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION DPHIO, TPIDEG, ONE
#endif
      DIMENSION  XC(*)
      REAL       X0(3)
      INTEGER    IDTYP(3,12)
      PARAMETER (TPIDEG=360,ONE=1)
      SAVE IDTYP
C.
      DATA  IDTYP / 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 2, 3, 1,
     +              2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 4, 3, 1, 1, 1,
     +              2, 3, 1, 2, 3, 1/
C.
C.    ------------------------------------------------------------------
C.
      JDIV  = LQ(JVO-1)
      ISH   = Q(JVO+2)
      IAXIS = Q(JDIV+1)
      IVOT  = Q(JDIV+2)
      JVOT  = LQ(JVOLUM-IVOT)
      IF (NLEVEL.LT.NLDEV(NLEVEL)) THEN
         JPAR = 0
      ELSE
*       (case with structure JVOLUM locally developed)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         IF (NLEVEL.EQ.NLDEV(NLEVEL)) GO TO 20
         DO 10 ILEV = NLDEV(NLEVEL), NLEVEL-1
            IF (IQ(JPAR+1).EQ.0) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
               IF (JPAR.EQ.0) GO TO 20
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
   10    CONTINUE
         GO TO 30
      ENDIF
*      (normal case)
   20 NDIV  = Q(JDIV+3)
      ORIG  = Q(JDIV+4)
      SDIV  = Q(JDIV+5)
*
   30 IDT = IDTYP(IAXIS,ISH)
      IF (IFL.NE.0) THEN
         IF (IDT.EQ.1) THEN
*
*         Division along X, Y or Z axis
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
*          Division along R axis
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
                  DO 40 IZ = 1, NZ-1
                     IF ((XC(3)-Q(IPT+3*IZ))*(XC(3)-Q(IPT+3*IZ+3))
     +               .LE.0.) THEN
                        IZSEC = IZ
                        IPT = IPT + 3 * IZSEC
                        GO TO 50
                     ENDIF
   40             CONTINUE
                  IN  = 0
                  GO TO 60
               ENDIF
   50          POR1 = (Q(IPT+3) - XC(3)) / (Q(IPT+3) - Q(IPT))
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
                     PH0 = PHI-Q(IPAR+1)
                     SG = SIGN(1.0,PH0)
                     PH0 = MOD( ABS(PH0), 360.0 )
                     IF(SG.LE.0.0) PH0 = 360.0-PH0
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
*          Division along Phi axis
*
            IF (XC(1).NE.0..OR.XC(2).NE.0.) THEN
               PHI = RADDEG * ATAN2 (XC(2), XC(1))
            ELSE
               PHI = 0.
            ENDIF
            DPHIO = PHI-ORIG
            SG = SIGN(ONE,DPHIO)
            DPHIO = MOD( ABS(DPHIO), TPIDEG)
            IF(SG.LE.0.0) DPHIO=TPIDEG-DPHIO
            IN = DPHIO/SDIV+1
         ELSE IF (IDT.EQ.4) THEN
*
*          Division along Theta axis
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
      ENDIF
*
   60 IF (IN.GT.NDIV) IN = 0
      IF (IN.LE.0) GO TO 999
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
      NL1    = NLEVEL
      NLEVEL = NLEVEL +1
      LVOLUM(NLEVEL) = IVOT
      NAMES(NLEVEL)  = IQ(JVOLUM+IVOT)
      NUMBER(NLEVEL) = IN
      LINDEX(NLEVEL) = IN
      LINMX(NLEVEL)  = NDIV
      GONLY(NLEVEL)  = GONLY(NL1)
      IF (LQ(LQ(JVOLUM-IVOT)).EQ.0) THEN
         NLDEV(NLEVEL) = NLDEV(NL1)
      ELSE
         NLDEV(NLEVEL) = NLEVEL
      ENDIF
*
      IF (IDT.EQ.1) THEN
         X0(1) = 0.0
         X0(2) = 0.0
         X0(3) = 0.0
         X0(IAXIS) = ORIG + (IN - 0.5) * SDIV
         IF (ISH.EQ.4.OR.(ISH.EQ.10.AND.IAXIS.NE.1)) THEN
            CALL GCENT (IAXIS, X0)
         ENDIF
         IF (IFL.NE.0) THEN
            XC(1) = XC(1) - X0(1)
            XC(2) = XC(2) - X0(2)
            XC(3) = XC(3) - X0(3)
         ENDIF
C*****  Code Expanded From Routine:  GTRMUL
C.
C.    ------------------------------------------------------------------
C.
         IF (GRMAT(10,NL1) .EQ. 0.0) THEN
            GTRAN(1,NLEVEL) = GTRAN(1,NL1) + X0(1)
            GTRAN(2,NLEVEL) = GTRAN(2,NL1) + X0(2)
            GTRAN(3,NLEVEL) = GTRAN(3,NL1) + X0(3)
            DO 70 I = 1, 10, 2
               GRMAT(I,NLEVEL) = GRMAT(I,NL1)
               GRMAT(I+1,NLEVEL) = GRMAT(I+1,NL1)
   70       CONTINUE
         ELSE
C
            DXTEM1 = X0(1)*GRMAT(1,NL1) + X0(2)*GRMAT(4,NL1) + X0(3)*
     +      GRMAT( 7,NL1)
            DXTEM2 = X0(1)*GRMAT(2,NL1) + X0(2)*GRMAT(5,NL1) + X0(3)*
     +      GRMAT( 8,NL1)
            DXTEM3 = X0(1)*GRMAT(3,NL1) + X0(2)*GRMAT(6,NL1) + X0(3)*
     +      GRMAT( 9,NL1)
            DO 80 I = 1, 10, 2
               GRMAT(I,NLEVEL) = GRMAT(I,NL1)
               GRMAT(I+1,NLEVEL) = GRMAT(I+1,NL1)
   80       CONTINUE
            GTRAN(1,NLEVEL) = GTRAN(1,NL1) + DXTEM1
            GTRAN(2,NLEVEL) = GTRAN(2,NL1) + DXTEM2
            GTRAN(3,NLEVEL) = GTRAN(3,NL1) + DXTEM3
         ENDIF
C*****  End of Code Expanded From Routine:  GTRMUL
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
         GTRAN(1,NLEVEL) = GTRAN(1,NL1)
         GRMAT(1,NLEVEL) = GRMAT(1,NL1)*CPHR + GRMAT(4,NL1)*SPHR
         GRMAT(4,NLEVEL) = GRMAT(4,NL1)*CPHR - GRMAT(1,NL1)*SPHR
         GRMAT(7,NLEVEL) = GRMAT(7,NL1)
         GTRAN(2,NLEVEL) = GTRAN(2,NL1)
         GRMAT(2,NLEVEL) = GRMAT(2,NL1)*CPHR + GRMAT(5,NL1)*SPHR
         GRMAT(5,NLEVEL) = GRMAT(5,NL1)*CPHR - GRMAT(2,NL1)*SPHR
         GRMAT(8,NLEVEL) = GRMAT(8,NL1)
         GTRAN(3,NLEVEL) = GTRAN(3,NL1)
         GRMAT(3,NLEVEL) = GRMAT(3,NL1)*CPHR + GRMAT(6,NL1)*SPHR
         GRMAT(6,NLEVEL) = GRMAT(6,NL1)*CPHR - GRMAT(3,NL1)*SPHR
         GRMAT(9,NLEVEL) = GRMAT(9,NL1)
         IF (IFL.NE.0) THEN
            XTT   = XC(1) * CPHR + XC(2) * SPHR
            XC(2) = XC(2) * CPHR - XC(1) * SPHR
            XC(1) = XTT
         ENDIF
         IF (PH0.EQ.0.0.AND.GRMAT(10,NL1).EQ.0.0) THEN
            GRMAT(10,NLEVEL) = 0.0
         ELSE
            GRMAT(10,NLEVEL) = 1.0
         ENDIF
         IF (ISH.EQ.11) IPSEC = 1
*
      ELSE
         GTRAN(1,NLEVEL) = GTRAN(1,NL1)
         GTRAN(2,NLEVEL) = GTRAN(2,NL1)
         GTRAN(3,NLEVEL) = GTRAN(3,NL1)
         DO 90 I = 1, 10, 2
            GRMAT(I,NLEVEL) = GRMAT(I,NL1)
            GRMAT(I+1,NLEVEL) = GRMAT(I+1,NL1)
   90    CONTINUE
      ENDIF
*
      IQ(JGPAR+NLEVEL) = NPAR
      LQ(JGPAR-NLEVEL) = JPAR
*                                                             END GMEDIV
  999 END

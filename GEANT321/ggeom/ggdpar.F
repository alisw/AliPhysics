*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GGDPAR (JVOM, IN, NVAR, LVAR, LVOM, NPAR, PAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GGDPAR (JVOM,IN,NVAR,LVAR,LVOM,NPAR*,PAR*)            *
C.    *                                                                *
C.    *   Computes the actual parameters for the IN'th division of the *
C.    *    mother volume at address JVOM                               *
C.    *   Returns the set of actual parameters in NPAR, PAR            *
C.    *                                                                *
C.    *   Called by : GGDVLP                                           *
C.    *   Author    : S.Banerjee                                       *
C.    *               (Original algorithms of A.McPherson)             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
C.
      PARAMETER (NPAMAX=50)
C.
      DIMENSION LVAR(*), PAR(*)
C.
      REAL      PARM(NPAMAX)
      INTEGER   LAX(10,12)
      SAVE LAX
C.
      DATA      LAX / 1, 2, 3, 7*0, 1, 1, 2, 3, 6*0, 1, 1, 2, 2, 3, 5*0,
     +                3, 0, 0, 2, 1, 1, 0, 2, 1, 1, 1, 1, 3, 7*0,
     +                1, 1, 3, 2, 2, 5*0, 3, 1, 1, 1, 1, 5*0,
     +                3, 1, 1, 1, 1, 2, 2, 3*0, 1, 1, 2, 2, 6*0,
     +                1, 2, 3, 7*0, 2, 2, 0, 0, 3, 1, 1, 3, 1, 1,
     +                2, 2, 0, 3, 1, 1, 3, 1, 1, 0/
C.
C.    ------------------------------------------------------------------
*
* *** Prepares parameters for mother, in PARM, for division, in PAR
*      and the division parameters
*
      JIN = LQ(JVOM-1)
      IF (LQ(JVOM).EQ.LVOM) THEN
*
*      Case when current volume is source of local development
*
         NPARM = Q(JVOM+5)
         CALL UCOPY (Q(JVOM+7), PARM, NPARM)
         NDIV  = Q(JIN+3)
         ORIG  = Q(JIN+4)
         STEP  = Q(JIN+5)
      ELSE
*
*      Other cases
*
         NPARM = IQ(LVOM+5)
         CALL UCOPY (Q(LVOM+6), PARM, NPARM)
         NDIV  = IQ(LVOM+1)
         ORIG  = Q(LVOM+2)
         STEP  = Q(LVOM+3)
      ENDIF
*
      ISHM  = Q(JVOM+2)
*
      ORI   = ORIG + (IN - 1) * STEP
      IAXIS = Q(JIN+1)
      IVO   = Q(JIN+2)
      JVO   = LQ(JVOLUM-IVO)
      ISH   = Q(JVO+2)
      NPAR  = Q(JVO+5)
*
* *** Prepare the division parameters
*
      IF (NPAR.GT.0) THEN
         CALL UCOPY (Q(JVO+7), PAR, NPAR)
      ELSE
         NPAR  = NPARM
         CALL UCOPY (PARM, PAR, NPARM)
*
*  **    Special treatment for phi divisions (when NPAR=0)
*
         IF ((ISHM.GE.5.AND.ISHM.LE.8.AND.IAXIS.EQ.2) .OR.
     +       (ISHM.GE.11.AND.ISHM.LE.12.AND.IAXIS.EQ.2) .OR.
     +       (ISHM.EQ.9.AND.IAXIS.EQ.3)) THEN
            IF (ISHM.EQ.5.OR.ISHM.EQ.7) THEN
               NPAR = NPARM + 2
               PAR(NPAR-1) = -0.5 * STEP
               PAR(NPAR)   =  0.5 * STEP
            ELSE IF (ISHM.EQ.6.OR.ISHM.EQ.8) THEN
               DP   = PARM(NPAR) - PARM(NPAR-1)
               IF (DP.LT.0.0) DP = DP + 360.0
               IF (ORIG-PARM(NPAR-1).LT.0.0) ORIG = ORIG + 360.0
               IF (ORIG-PARM(NPAR-1).GT.DP ) GO TO 910
               DP   = PARM(NPAR) - ORIG
               IF (DP.LT.0.0) DP = DP + 360.0
               PAR(NPAR-1) = -0.5 * DP / NDIV
               PAR(NPAR  ) =  0.5 * DP / NDIV
            ELSE IF (ISHM.EQ.11.OR.ISHM.EQ.12) THEN
               IF (ISHM.EQ.11) NDIV = PARM(3)
               STEP = PARM(2) / NDIV
               PAR(1) = -0.5 * STEP
               PAR(2) = STEP
               PAR(3) = 1.
            ELSE IF (ISHM.EQ.9) THEN
               DP   = PARM(6) - PARM(5)
               IF (DP.LT.0.0) DP = DP + 360.0
               IF (ORIG-PARM(5).LT.0.0) ORIG = ORIG + 360.0
               IF (ORIG-PARM(5).GT.DP ) GO TO 910
               DP   = PARM(6) - ORIG
               IF (DP.LT.0.0) DP = DP + 360.0
               PAR(5) = -0.5 * DP / NDIV
               PAR(6) =  0.5 * DP / NDIV
            ENDIF
         ENDIF
      ENDIF
      IF (NVAR.LE.0)  GO TO 999
*
* *** Compute the actual parameters
*
      IF (ISHM.EQ.1) THEN
*
* BOX
*
         IF (ISH.EQ.1) THEN
            DO 10 I = 1, NVAR
               IAX = LVAR(I)
               IF (IAX.EQ.IAXIS) THEN
                  PAR(IAX) = 0.5 *STEP
               ELSE
                  PAR(IAX) = PARM(IAX)
               ENDIF
   10       CONTINUE
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.2) THEN
*
* TRD1
*
         IF (ISH.EQ.2) THEN
            DO 20 I = 1, NVAR
               IAX = LVAR(I)
               IF (LAX(IAX,ISH).EQ.IAXIS) THEN
                  PAR(IAX) = 0.5 * STEP
               ELSE IF (LAX(IAX,ISH).EQ.1.AND.IAXIS.EQ.3) THEN
                  ZZ   = ORI + PARM(4)
                  DXDZ = 0.5 * (PARM(2) - PARM(1)) / PARM(4)
                  IF (IAX.EQ.2) ZZ = ZZ + STEP
                  PAR(IAX) = PARM(1) + DXDZ * ZZ
               ELSE
                  PAR(IAX) = PARM(IAX)
               ENDIF
   20       CONTINUE
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.3) THEN
*
* TRD2
*
         IF (ISH.EQ.3.AND.IAXIS.EQ.3) THEN
            DO 30 I = 1, NVAR
               IAX = LVAR(I)
               IF (LAX(IAX,ISH).EQ.IAXIS) THEN
                  PAR(IAX) = 0.5 * STEP
               ELSE
                  IP1  = 2 * LAX(IAX,ISH) - 1
                  IP2  = IP1 + 1
                  ZZ   = ORI + PARM(5)
                  DXDZ = 0.5 * (PARM(IP2) - PARM(IP1)) / PARM(5)
                  IF (IAX.EQ.IP2) ZZ = ZZ + STEP
                  PAR(IAX) = PARM(IP1) + DXDZ * ZZ
               ENDIF
   30       CONTINUE
         ELSE IF (ISH.EQ.3) THEN
            GO TO 910
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.4) THEN
*
* TRAP
*
         IF (ISH.EQ.4.AND.IAXIS.NE.1) THEN
         IF(IAXIS.EQ.3) THEN
            DO 40 I = 1, NVAR
               IAX = LVAR(I)
               IF (IAX.EQ.1) THEN
                  PAR(IAX) = 0.5 * STEP
               ELSE IF (IAX.LE.6) THEN
                  ZZ   = ORI + PARM(1)
                  DPDZ = 0.5 * (PARM(IAX+4) - PARM(IAX)) / PARM(1)
                  PAR(IAX) = PARM(IAX) + DPDZ * ZZ
               ELSE
                  ZZ   = ORI + PARM(1) + STEP
                  DPDZ = 0.5 * (PARM(IAX) - PARM(IAX-4)) / PARM(1)
                  PAR(IAX) = PARM(IAX-4) + DPDZ * ZZ
               ENDIF
   40       CONTINUE
            HTAL = PARM(8) * PARM(11)
            HTAH = PARM(4) * PARM(7)
            ZZ1  = 0.5 * (ORI + PARM(1)) / PARM(1)
            ZZ2  = 0.5 * (ORI + PARM(1) + STEP) / PARM(1)
            PAR(7) = (HTAL*(1.0-ZZ1) + HTAH*ZZ1) / PARM(4)
            PAR(11)= (HTAL*(1.0-ZZ2) + HTAH*ZZ2) / PARM(8)
            CALL GNOTR1(PAR)
            ELSE IF (IAXIS.EQ.2) THEN
               DO 41 I = 1, NVAR
               PAR(LVAR(I)) = PARM(LVAR(I))
  41           CONTINUE
               DXDY1 = 0.5*(PARM(6)-PARM(5))/PARM(4)
               DXDY2 = 0.5*(PARM(10)-PARM(9))/PARM(8)
               DXM1  = 0.5*(PARM(6)+PARM(5))
               DXM2  = 0.5*(PARM(10)+PARM(9))
               DXH1  = DXM1+(ORI+STEP)*DXDY1
               DXH2  = DXM2+(ORI+STEP)*DXDY2
               DXL1  = DXM1+ORI*DXDY1
               DXL2  = DXM2+ORI*DXDY2
               PAR(5) = DXL1
               PAR(6) = DXH1
               PAR(9) = DXL2
               PAR(10) = DXH2
            ENDIF
         ELSE IF (ISH.EQ.4) THEN
            GO TO 910
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.5 .OR. ISHM.EQ.6) THEN
*
* TUBE or TUBS
*
         IF (ISH.EQ.5 .OR. ISH.EQ.6) THEN
            DO 50 I = 1, NVAR
               IAX = LVAR(I)
               IF (LAX(IAX,ISH).EQ.IAXIS) THEN
                  IF (IAXIS.EQ.3) THEN
                     PAR(IAX) = 0.5 * STEP
                  ELSE IF (IAXIS.EQ.1) THEN
                     IF (IAX.EQ.1) THEN
                        PAR(IAX) = ORI
                     ELSE
                        PAR(IAX) = ORI + STEP
                     ENDIF
                  ELSE
                     GO TO 910
                  ENDIF
               ELSE
                  PAR(IAX) = PARM(IAX)
               ENDIF
   50       CONTINUE
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.7 .OR. ISHM.EQ.8) THEN
*
* CONE or CONS
*
         IF (ISH.EQ.7 .OR. ISH.EQ.8) THEN
            DO 60 I = 1, NVAR
               IAX = LVAR(I)
               IF (LAX(IAX,ISH).EQ.IAXIS .AND. IAXIS.EQ.3) THEN
                  PAR(IAX) = 0.5 * STEP
               ELSE IF (IAXIS.EQ.3.AND.IAX.GT.1.AND.IAX.LT.6) THEN
                  IF (IAX.EQ.2.OR.IAX.EQ.4) THEN
                     DP   = 0.5 * (PARM(4) - PARM(2)) / PARM(1)
                     PM   = 0.5 * (PARM(4) + PARM(2))
                  ELSE
                     DP   = 0.5 * (PARM(5) - PARM(3)) / PARM(1)
                     PM   = 0.5 * (PARM(5) + PARM(3))
                  ENDIF
                  IF (IAX.GT.3) THEN
                     DZ   = ORI + STEP
                  ELSE
                     DZ   = ORI
                  ENDIF
                  PAR(IAX) = PM + DP * DZ
               ELSE IF (IAXIS.EQ.1.AND.LAX(IAX,ISH).EQ.IAXIS) THEN
                  IF (IAX.EQ.2) THEN
                     PAR(IAX) = ORI
                  ELSE IF (IAX.EQ.3) THEN
                     PAR(IAX) = ORI + STEP
                  ELSE IF (IAX.EQ.4) THEN
                     PAR(IAX) = ORI * PARM(IAX) / PARM(2)
                  ELSE
                     PAR(IAX) = (ORI + STEP) * PARM(IAX) / PARM(3)
                  ENDIF
               ELSE
                  PAR(IAX) = PARM(IAX)
               ENDIF
   60       CONTINUE
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.9) THEN
*
* SPHE
*
         IF (ISH.EQ.9) THEN
            DO 70 I = 1, NVAR
               IAX = LVAR(I)
               IF (LAX(IAX,ISH).EQ.IAXIS) THEN
                  IF (MOD(IAX,2).NE.0) THEN
                     PAR(IAX) = ORI
                  ELSE
                     PAR(IAX) = ORI + STEP
                  ENDIF
               ELSE
                  PAR(IAX) = PARM(IAX)
               ENDIF
   70       CONTINUE
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.10) THEN
*
* PARA
*
         IF (ISH.EQ.10) THEN
            DO 80 I = 1, NVAR
               IAX = LVAR(I)
               IF (LAX(IAX,ISH).EQ.IAXIS) THEN
                  PAR(IAX) = 0.5 * STEP
               ENDIF
   80       CONTINUE
         ELSE
            GO TO 900
         ENDIF
*
      ELSE IF (ISHM.EQ.11 .OR. ISHM.EQ.12) THEN
*
* PGON or PCON
*
         IF (ISH.EQ.ISHM) THEN
            IF (ISH.EQ.11) THEN
               IPNZ = 4
            ELSE
               IPNZ = 3
            ENDIF
            NZ   = PAR(IPNZ)
            NZ1  = PARM(IPNZ)
            IF (NZ.EQ.NZ1) THEN
               IF (IAXIS.EQ.1) THEN
                  RMN  = ORI
                  RMX  = ORI + STEP
                  DO 90 I = 1, NZ
                     IAX1 = IPNZ + 3*I - 1
                     IAX2 = IAX1 + 1
                     IF (I.EQ.1) THEN
                        PAR(IAX1) = RMN
                        PAR(IAX2) = RMX
                     ELSE
                        PAR(IAX1) = RMN * PARM(IAX1) / PARM(IPNZ+2)
                        PAR(IAX2) = RMX * PARM(IAX2) / PARM(IPNZ+3)
                     ENDIF
   90             CONTINUE
               ELSE IF (IAXIS.EQ.2) THEN
                  DO 100 I = 1, NVAR
                     IAX = LVAR(I)
                     PAR(IAX) = PARM(IAX)
  100             CONTINUE
               ELSE IF (NZ.EQ.2.AND.IAXIS.EQ.3) THEN
                  ZL   = PARM(IPNZ+1)
                  ZH   = PARM(IPNZ+4)
                  DZ   = ZH - ZL
                  DRMIDZ = (PARM(IPNZ+5)-PARM(IPNZ+2))/DZ
                  DRMADZ = (PARM(IPNZ+6)-PARM(IPNZ+3))/DZ
                  PAR(IPNZ+1) = -0.5 * STEP
                  PAR(IPNZ+4) =  0.5 * STEP
                  DO 110 I = 1, NVAR
                     IAX = LVAR(I)-IPNZ
                     IF(IAX.EQ.2) THEN
                        RAD = PARM(IPNZ+2)+(IN-1)*STEP*DRMIDZ
                     ELSEIF (IAX.EQ.3) THEN
                        RAD = PARM(IPNZ+3)+(IN-1)*STEP*DRMADZ
                     ELSEIF (IAX.EQ.5) THEN
                        RAD = PARM(IPNZ+2)+IN*STEP*DRMIDZ
                     ELSEIF (IAX.EQ.6) THEN
                        RAD = PARM(IPNZ+3)+IN*STEP*DRMADZ
                     ENDIF
                     PAR(IPNZ+IAX) = RAD
  110             CONTINUE
               ELSE
                  GO TO 920
               ENDIF
            ELSE
               GO TO 920
            ENDIF
         ELSE
            GO TO 900
         ENDIF
*
      ELSE
         GO TO 900
      ENDIF
*
      GO TO 999
*
  900 WRITE (CHMAIL, 1001) ISH, ISHM
      GO TO 990
*
  910 WRITE (CHMAIL, 1002) ISH, ISHM, IAXIS
      GO TO 990
*
  920 WRITE (CHMAIL, 1003) ISH, NZ, ISHM, NZ1
*
  990 CALL GMAIL( 0, 0)
      IEORUN = 1
*
 1001 FORMAT (' GGDPAR : Not accepted ISH,ISHM=',2I5)
 1002 FORMAT (' GGDPAR : Not accepted ISH,ISHM,IAXIS=',3I5)
 1003 FORMAT (' GGDPAR : Not accepted ISH,NZ,ISHM,NZ1=',4I5)
*                                                             END GGDPAR
  999 END

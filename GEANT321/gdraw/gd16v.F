*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_CG)
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GD16V(IWOFFS,SHADE)
************************************************************************
*                                                                      *
*     Name: GD16V                                                      *
*     Author: S.Giani                            Date:    18.07.91     *
*                                                Revised: 1992         *
*                                                                      *
*     Function: Vizualisation of WIRE-object                           *
*               Search and construction of faces' visible portions     *
*               Surface filling for WIRE-object                        *
*               Two scan line algorithms                               *
*               Surface shading                                        *
*                                                                      *
*     References: CGVHED                                               *
*                                                                      *
*     Input: Q(JCG+IWOFFS+*) - WIRE-object                             *
*               Q(JCG+*) - HIDE-structure                              *
*                                                                      *
*     Output: none                                                     *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/cggpar.inc"
#include "geant321/cghpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/cgctra.inc"
#include "geant321/cgcedg.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcspee.inc"
#include "geant321/gconst.inc"
#include "geant321/gchiln.inc"
*
      INTEGER SHADE(*)
      REAL P1(3),P2(3),AX(2),AY(2)
*SG
      DIMENSION ISCFAC(500),SCXYZ(3,500),COSDIR(4)
      DIMENSION CX(500),CY(500),CZ(500)
      DIMENSION TESTX(2),TESTY(2)
      DIMENSION VVX(500),VVY(500)
*
      DIMENSION AROT(4,4),ZROT(4,4),CXROT(500),CYROT(500),CZROT(500)
      DIMENSION PXINT(50),PZINT(50),XPINT(50),YPINT(50),ZPINT(50),RRR(4)
*SG
      LLEP=ABS(LEP)
*****                                 CALL IGSET('SYNC',1.)
      IF (Q(JCG+IWOFFS+KCGSIZ) .LE. 0.)   THEN
         WRITE(CHMAIL,10000)
         GOTO 400
      ENDIF
      IF (Q(JCG+KHSIZE) .LE. LHHEAD)  THEN
         WRITE(CHMAIL,10100)
         GOTO 400
      ENDIF
*SG
      LINBUF=0
      DO 10  IJ=1,2
         TESTX(IJ)=20.
         TESTY(IJ)=20.
   10 CONTINUE
      LINFIL=IBITS(LINATT,13,3)
      IF(LINFIL.NE.0)CALL GDRAWV(TESTX,TESTY,-1)
      IJKLMN=0
      IACTUN=0
      LPARZ=0
      LFULL=0
      IN=0
      KKK=0
      LMN=0
      NSI=0
      IGEN=0
*SG
      NT = Q(JCG+KHNT)
      IF (NT.LE.0 .OR. NT.GT.NTMAX) GOTO 400
      NEDGE  = Q(JCG+IWOFFS+KCGNF)
*      IF(NEDGE.EQ.0)GOTO 49
      JXYZ1   = Q(JCG+KHJXYZ)
      JFA    = Q(JCG+KHJFA)
      JPFA   = Q(JCG+KHJPFA)
      JDFA   = Q(JCG+KHJDFA)
      JTRE   = Q(JCG+KHJTRE)
      JSTA   = Q(JCG+KHJSTA)
      JALE   = Q(JCG+KHJALE)
      JARI   = Q(JCG+KHJARI)
      NFACE  = Q(JCG+KHNFAC)
      DO 230 NE=1,NEDGE
         IF(LINFIL.NE.0)THEN
            IF(LLEP.NE.1)THEN
               DO 180 MMM=2,SHADE(1)
                  IF(NE.EQ.(SHADE(MMM)+1).AND.NE.NE.1)THEN
***                     IF((LFULL+LPARZ).EQ.0)GOTO 260
                     CALL CGFAC2(CX,CY,CZ,LMN)
                     ISCFAC(1)=LMN*.5
                     DO 20  I=2,LMN+1
                        ISCFAC(I)=I-1
                        SCXYZ(1,I-1)=CX(I-1)
                        SCXYZ(2,I-1)=CY(I-1)
                        SCXYZ(3,I-1)=CZ(I-1)
   20                CONTINUE
                     CALL CGHPLA(ISCFAC,SCXYZ,COSDIR)
                     XCOSX=(SIN(GTHETA*DEGRAD))*(COS(GPHI*DEGRAD))
                     YCOSY=(SIN(GTHETA*DEGRAD))*(SIN(GPHI*DEGRAD))
                     ZCOSZ=COS(GTHETA*DEGRAD)
                     PROSCA=(COSDIR(1)*XCOSX)+(COSDIR(2)*YCOSY)+
     +               (COSDIR(3)*ZCOSZ)
                     PROSCB=(2.*(PROSCA**2))-1.
                     APROSC=PROSCB
*            IF(APROSC.LT.-1.)PRINT *,'GT1'
                     IF(LFULL.EQ.(IJKLMN/2).AND.
     +       (LINFIL.EQ.2.OR.LINFIL.EQ.3))THEN
                        CALL CGFACO(VVX,VVY,IJKLMN,LFULL,APROSC)
                     ELSE
                        IF(ISWIT(10).EQ.100)THEN
                           PNX=COSDIR(1)
                           PNY=COSDIR(2)
                           PNZ=COSDIR(3)
                           APNZ=ABS(1.-PNZ)
                           IF(APNZ.LT..0001)THEN
                              AROT(1,1)=1.
                              AROT(1,2)=0.
                              AROT(1,3)=0.
                              AROT(2,1)=0.
                              AROT(2,2)=1.
                           ELSE
                              AROT(1,1)=PNX*PNZ/SQRT(PNX**2+PNY**2)
                              AROT(1,2)=PNY*PNZ/SQRT(PNX**2+PNY**2)
                              AROT(1,3)=-SQRT(PNX**2+PNY**2)
                              AROT(2,1)=-PNY/SQRT(PNX**2+PNY**2)
                              AROT(2,2)=PNX/SQRT(PNX**2+PNY**2)
                           ENDIF
                           AROT(2,3)=0.
                           AROT(3,1)=PNX
                           AROT(3,2)=PNY
                           AROT(3,3)=PNZ
                        ELSE
                           DO 40  IHH=1,3
                              DO 30  JHH=1,4
                                 AROT(IHH,JHH)=TSCRN(JHH,IHH,NT)
   30                         CONTINUE
   40                      CONTINUE
                           AROT(4,1)=0.
                           AROT(4,2)=0.
                           AROT(4,3)=0.
                           AROT(4,4)=1.
                        ENDIF
                        IF(ISWIT(10).EQ.100)THEN
                           ZROT(1,1)=AROT(1,1)
                           ZROT(2,1)=AROT(1,2)
                           ZROT(3,1)=AROT(1,3)
                           ZROT(1,2)=AROT(2,1)
                           ZROT(2,2)=AROT(2,2)
                           ZROT(3,2)=AROT(2,3)
                           ZROT(1,3)=AROT(3,1)
                           ZROT(2,3)=AROT(3,2)
                           ZROT(3,3)=AROT(3,3)
                        ELSE
                           DO 60  IHH=1,4
                              DO 50  JHH=1,4
                                 ZROT(IHH,JHH)=AROT(IHH,JHH)
   50                         CONTINUE
   60                      CONTINUE
                           CALL RINV(4,ZROT,4,RRR,IFAIL)
                        ENDIF
                        YROTMI=100000.
                        YROTMA=-100000.
                        DO 70 M=1,LMN
                           IF(ISWIT(10).EQ.100)THEN
                              CXROT(M)=AROT(1,1)*CX(M)+ AROT(1,2)*CY(M)
     +                        +AROT(1,3)*CZ(M)
                              CYROT(M)=AROT(2,1)*CX(M)+ AROT(2,2)*CY(M)
     +                        +AROT(2,3)*CZ(M)
                              CZROT(M)=AROT(3,1)*CX(M)+ AROT(3,2)*CY(M)
     +                        +AROT(3,3)*CZ(M)
                           ELSE
                              CXROT(M)=AROT(1,4)+AROT(1,1)*CX(M)+
     +                        AROT(1,2)*CY(M)+AROT(1,3)*CZ(M)
                              CYROT(M)=AROT(2,4)+AROT(2,1)*CX(M)+
     +                        AROT(2,2)*CY(M)+AROT(2,3)*CZ(M)
                              CZROT(M)=AROT(3,4)+AROT(3,1)*CX(M)+
     +                        AROT(3,2)*CY(M)+AROT(3,3)*CZ(M)
                           ENDIF
                           IF(CYROT(M).LT.YROTMI)YROTMI=CYROT(M)
                           IF(CYROT(M).GT.YROTMA)YROTMA=CYROT(M)
   70                   CONTINUE
                        IF(LINFIL.EQ.1)THEN
                           RINULI=.01
                        ELSEIF(LINFIL.EQ.2)THEN
                           RINULI=.1
                        ELSEIF(LINFIL.EQ.3)THEN
                           RINULI=.05
                        ELSEIF(LINFIL.EQ.4)THEN
                           RINULI=.05
                        ELSEIF(LINFIL.EQ.5)THEN
                           RINULI=.02
                        ELSEIF(LINFIL.EQ.6)THEN
                           RINULI=.01
                        ELSEIF(LINFIL.EQ.7)THEN
                           RINULI=.005
                        ENDIF
                        NYROTM=(YROTMA-YROTMI)/RINULI
                        YROTST=YROTMI
                        DO 170 MM=1,NYROTM-1
                           YROTST=YROTST+RINULI
                           JFK=0
                           DO 80 MMI=1,LMN-1,2
                              IF(CXROT(MMI).EQ.CXROT(MMI+1))THEN
                                 IF(MIN(CYROT(MMI),CYROT(MMI+1)).LT.
     +                           YROTST.AND.MAX(CYROT(MMI),CYROT(MMI+1)
     +                           ) .GT.YROTST)THEN
                                    JFK=JFK+1
                                    PXINT(JFK)=CXROT(MMI)
                                    IF(CZROT(MMI).EQ.CZROT(MMI+1))THEN
                                    PZINT(JFK)=CZROT(MMI)
                                 ELSEIF(CYROT(MMI).NE.CYROT(MMI+1))
     +                           THEN
                                    AAAZ=(CYROT(MMI)-CYROT(MMI+1))/
     +                              (CZROT(MMI)-CZROT(MMI+1))
                                    BBBZ=(CZROT(MMI)*CYROT(MMI+1)-
     +                              CYROT(MMI)*CZROT(MMI+1))/ (CZROT(MM
     +                              I)-CZROT(MMI+1))
                                    PZINT(JFK)=(YROTST-BBBZ)/AAAZ
                                    ENDIF
                                 ENDIF
                              ELSEIF(CYROT(MMI).NE.CYROT(MMI+1))THEN
                                 IF(MIN(CYROT(MMI),CYROT(MMI+1)).LT.
     +                           YROTST.AND.MAX(CYROT(MMI),CYROT(MMI+1)
     +                           ) .GT.YROTST)THEN
                                    AAA=(CYROT(MMI)-CYROT(MMI+1))/
     +                              (CXROT(MMI)-CXROT(MMI+1))
                                    BBB=(CXROT(MMI)*CYROT(MMI+1)-
     +                              CYROT(MMI)*CXROT(MMI+1))/ (CXROT(MM
     +                              I)-CXROT(MMI+1))
                                    JFK=JFK+1
                                    PXINT(JFK)=(YROTST-BBB)/AAA
                                    IF(CZROT(MMI).EQ.CZROT(MMI+1))THEN
                                    PZINT(JFK)=CZROT(MMI)
                                 ELSEIF(CYROT(MMI).NE.CYROT(MMI+1))
     +                           THEN
                                    AAAZ=(CYROT(MMI)-CYROT(MMI+1))/
     +                              (CZROT(MMI)-CZROT(MMI+1))
                                    BBBZ=(CZROT(MMI)*CYROT(MMI+1)-
     +                              CYROT(MMI)*CZROT(MMI+1))/ (CZROT(MM
     +                              I)-CZROT(MMI+1))
                                    PZINT(JFK)=(YROTST-BBBZ)/AAAZ
                                    ENDIF
                                 ENDIF
                              ENDIF
   80                      CONTINUE
                           DO 100 JM=1,JFK-1
                              DO 90 KM=JM+1,JFK
                                 IF(PXINT(JM).LT.PXINT(KM))THEN
                                    TMPM=PXINT(JM)
                                    TMPZ=PZINT(JM)
                                    PXINT(JM)=PXINT(KM)
                                    PZINT(JM)=PZINT(KM)
                                    PXINT(KM)=TMPM
                                    PZINT(KM)=TMPZ
                                 ENDIF
   90                         CONTINUE
  100                      CONTINUE
                           RJFK=JFK*.5
                           IRJFK=RJFK
*                          IF((RJFK-IRJFK).GT..1)PRINT *,'Odd !'
                           DO 110 MR=1,JFK
                              IF(ISWIT(10).EQ.100)THEN
                                 XPINT(MR)=ZROT(1,1)*PXINT(MR)+ ZROT(1,
     +                           2)*YROTST+ZROT(1,3)*CZROT(1)
                                 YPINT(MR)=ZROT(2,1)*PXINT(MR)+ ZROT(2,
     +                           2)*YROTST+ZROT(2,3)*CZROT(1)
                                 ZPINT(MR)=ZROT(3,1)*PXINT(MR)+ ZROT(3,
     +                           2)*YROTST+ZROT(3,3)*CZROT(1)
                              ELSE
                                 XPINT(MR)=ZROT(1,4)+ZROT(1,1)*
     +                           PXINT(MR)+ ZROT(1,2)*YROTST+ZROT(1,3)*
     +                           PZINT(MR)
                                 YPINT(MR)=ZROT(2,4)+ZROT(2,1)*
     +                           PXINT(MR)+ ZROT(2,2)*YROTST+ZROT(2,3)*
     +                           PZINT(MR)
                                 ZPINT(MR)=ZROT(3,4)+ZROT(3,1)*
     +                           PXINT(MR)+ ZROT(3,2)*YROTST+ZROT(3,3)*
     +                           PZINT(MR)
                              ENDIF
  110                      CONTINUE
                           IF(LINFIL.GT.1)THEN
                              LINCOL=IBITS(LINATT,16,8)
                              CALL GDSHAD(LINCOL,APROSC)
                           ENDIF
                           DO 160 MZ=1,JFK,2
                              P1(1)=XPINT(MZ)
                              P1(2)=YPINT(MZ)
                              P1(3)=ZPINT(MZ)
                              P2(1)=XPINT(MZ+1)
                              P2(2)=YPINT(MZ+1)
                              P2(3)=ZPINT(MZ+1)
                              CALL CGVEDG(NT,P1,P2,IVIS)
                              IF (IVIS .LT. 0) GOTO 150
                              IF (NFACE .EQ. 0) GOTO 140
                              CALL CGVHED( Q(JCG+JXYZ1),IQ(JCG+JFA),
     +                        IQ(JCG+JPFA), Q(JCG+ JDFA),IQ(JCG+JTRE),
     +                        IQ(JCG+JSTA), Q(JCG+JALE),Q(JCG+JARI))
*           D R A W   E D G E
                              IF (NPART) 150,140,120
  120                         DO 130 I=1,NPART+1
                                 IF (I .EQ. 1) T1 = 0.
                                 IF (I .NE. 1) T1 = TEND(I-1)
                                 IF (I .NE. NPART+1) T2 = TSTRT(I)
                                 IF (I .EQ. NPART+1) T2 = 1.
                                 IF (T2-T1 .LT. TDEL) GOTO 130
                                 AX(1) = XA + T1*XDELT
                                 AY(1) = YA + T1*YDELT
                                 AX(2) = XA + T2*XDELT
                                 AY(2) = YA + T2*YDELT
                                 IF(LINFIL.GT.1.AND.IDVIEW.NE.0)THEN
                                    LLLINA=LINATT
                                    CALL MVBITS(LINCOL,0,8,LINATT,16)
                                    CALL GDRAWV(AX,AY,2)
                                    LINATT=LLLINA
                                 ELSE
                                    CALL GDRAWV(AX,AY,2)
                                 ENDIF
  130                         CONTINUE
                              GOTO 150
*           L I N E   I S   F U L L Y   V I S I B L E
  140                         CONTINUE
                              AX(1) = AA(1)
                              AY(1) = AA(2)
                              AX(2) = BB(1)
                              AY(2) = BB(2)
                              IF(LINFIL.GT.1.AND.IDVIEW.NE.0)THEN
                                 LLLINA=LINATT
                                 CALL MVBITS(LINCOL,0,8,LINATT,16)
                                 CALL GDRAWV(AX,AY,2)
                                 LINATT=LLLINA
                              ELSE
                                 CALL GDRAWV(AX,AY,2)
                              ENDIF
  150                         CONTINUE
  160                      CONTINUE
  170                   CONTINUE
                     ENDIF
*  Resetting counters for next face
                     IJKLMN=0
                     LPARZ=0
                     LFULL=0
                     IN=0
                     KKK=0
                     LMN=0
                     NSI=0
                     IGEN=0
                     GOTO 190
                  ENDIF
  180          CONTINUE
            ENDIF
         ENDIF
  190    CONTINUE
         J = LCGHEA + (NE-1)*LCGEDG
*         IEA = Q(JCG+IWOFFS+J+KCGAE)
         P1(1) = Q(JCG+IWOFFS+J+KCGX1)
         P1(2) = Q(JCG+IWOFFS+J+KCGY1)
         P1(3) = Q(JCG+IWOFFS+J+KCGZ1)
         P2(1) = Q(JCG+IWOFFS+J+KCGX2)
         P2(2) = Q(JCG+IWOFFS+J+KCGY2)
         P2(3) = Q(JCG+IWOFFS+J+KCGZ2)
         IF(LINFIL.NE.0)THEN
            LMN=LMN+1
*          IF(LMN.GE.500)THEN
*            PRINT *,LMN,'=LMN'
*          ENDIF
            CX(LMN)=P1(1)
            CY(LMN)=P1(2)
            CZ(LMN)=P1(3)
            LMN=LMN+1
*          IF(LMN.GE.500)THEN
*            PRINT *,LMN,'=LMN'
*          ENDIF
            CX(LMN)=P2(1)
            CY(LMN)=P2(2)
            CZ(LMN)=P2(3)
         ENDIF
         CALL CGVEDG(NT,P1,P2,IVIS)
         IF(LINFIL.NE.0)THEN
            IJKLMN=IJKLMN+1
*              IF(IJKLMN.GE.500)THEN
*                PRINT *,IJKLMN,'=IJKLMN'
*              ENDIF
            VVX(IJKLMN)=AA(1)
            VVY(IJKLMN)=AA(2)
            IJKLMN=IJKLMN+1
*              IF(IJKLMN.GE.500)THEN
*                PRINT *,IJKLMN,'=IJKLMN'
*              ENDIF
            VVX(IJKLMN)=BB(1)
            VVY(IJKLMN)=BB(2)
         ENDIF
         IF (IVIS .LT. 0) GOTO 230
         IF (NFACE .EQ. 0) GOTO 220
         CALL CGVHED( Q(JCG+JXYZ1),IQ(JCG+JFA),IQ(JCG+JPFA), Q(JCG+
     +   JDFA),IQ(JCG+JTRE),IQ(JCG+JSTA), Q(JCG+JALE),Q(JCG+JARI))
*           D R A W   E D G E
         IF (NPART) 230,220,200
  200    DO 210 I=1,NPART+1
            IF (I .EQ. 1) T1 = 0.
            IF (I .NE. 1) T1 = TEND(I-1)
            IF (I .NE. NPART+1) T2 = TSTRT(I)
            IF (I .EQ. NPART+1) T2 = 1.
            IF (T2-T1 .LT. TDEL) GOTO 210
            AX(1) = XA + T1*XDELT
            AY(1) = YA + T1*YDELT
            AX(2) = XA + T2*XDELT
            AY(2) = YA + T2*YDELT
            IF(LINFIL.NE.0)THEN
               LTY=IBITS(LINATT,10,3)
               CALL MVBITS(6,0,3,LINATT,10)
***               IOLEP=LEP
***               LEP=11
               CALL GDRAWV(AX,AY,2)
***               LEP=IOLEP
               CALL MVBITS(LTY,0,3,LINATT,10)
            ELSE
               CALL GDRAWV(AX,AY,2)
            ENDIF
  210    CONTINUE
         LPARZ=LPARZ+1
         GOTO 230
*           L I N E   I S   F U L L Y   V I S I B L E
  220    CONTINUE
         AX(1) = AA(1)
         AY(1) = AA(2)
         AX(2) = BB(1)
         AY(2) = BB(2)
         IF(LINFIL.NE.0)THEN
            LTY=IBITS(LINATT,10,3)
            CALL MVBITS(6,0,3,LINATT,10)
***            IOLEP=LEP
***            LEP=11
            CALL GDRAWV(AX,AY,2)
***            LEP=IOLEP
            CALL MVBITS(LTY,0,3,LINATT,10)
         ELSE
            CALL GDRAWV(AX,AY,2)
         ENDIF
         LFULL=LFULL+1
  230 CONTINUE
*SG
      IF(LINFIL.NE.0)THEN
***         IF((LFULL+LPARZ).EQ.0)GOTO 555
         CALL CGFAC2(CX,CY,CZ,LMN)
         ISCFAC(1)=LMN*.5
         DO 240  I=2,LMN+1
            ISCFAC(I)=I-1
            SCXYZ(1,I-1)=CX(I-1)
            SCXYZ(2,I-1)=CY(I-1)
            SCXYZ(3,I-1)=CZ(I-1)
  240    CONTINUE
         CALL CGHPLA(ISCFAC,SCXYZ,COSDIR)
         XCOSX=(SIN(GTHETA*DEGRAD))*(COS(GPHI*DEGRAD))
         YCOSY=(SIN(GTHETA*DEGRAD))*(SIN(GPHI*DEGRAD))
         ZCOSZ=COS(GTHETA*DEGRAD)
         PROSCA=(COSDIR(1)*XCOSX)+(COSDIR(2)*YCOSY)+ (COSDIR(3)*ZCOSZ)
         PROSCB=(2.*(PROSCA**2))-1.
         APROSC=PROSCB
*             IF(APROSC.GT.1)PRINT *,'GT1'
         IF(LFULL.EQ.IJKLMN/2.AND.
     +       (LINFIL.EQ.2.OR.LINFIL.EQ.3))THEN
            CALL CGFACO(VVX,VVY,IJKLMN,LFULL,APROSC)
         ELSE
            IF(ISWIT(10).EQ.100)THEN
               PNX=COSDIR(1)
               PNY=COSDIR(2)
               PNZ=COSDIR(3)
               APNZ=ABS(1.-PNZ)
               IF(APNZ.LT..0001)THEN
                  AROT(1,1)=1.
                  AROT(1,2)=0.
                  AROT(1,3)=0.
                  AROT(2,1)=0.
                  AROT(2,2)=1.
               ELSE
                  AROT(1,1)=PNX*PNZ/SQRT(PNX**2+PNY**2)
                  AROT(1,2)=PNY*PNZ/SQRT(PNX**2+PNY**2)
                  AROT(1,3)=-SQRT(PNX**2+PNY**2)
                  AROT(2,1)=-PNY/SQRT(PNX**2+PNY**2)
                  AROT(2,2)=PNX/SQRT(PNX**2+PNY**2)
               ENDIF
               AROT(2,3)=0.
               AROT(3,1)=PNX
               AROT(3,2)=PNY
               AROT(3,3)=PNZ
            ELSE
               DO 260 IHH=1,3
                  DO 250 JHH=1,4
                     AROT(IHH,JHH)=TSCRN(JHH,IHH,NT)
  250             CONTINUE
  260          CONTINUE
               AROT(4,1)=0.
               AROT(4,2)=0.
               AROT(4,3)=0.
               AROT(4,4)=1.
            ENDIF
            IF(ISWIT(10).EQ.100)THEN
               ZROT(1,1)=AROT(1,1)
               ZROT(2,1)=AROT(1,2)
               ZROT(3,1)=AROT(1,3)
               ZROT(1,2)=AROT(2,1)
               ZROT(2,2)=AROT(2,2)
               ZROT(3,2)=AROT(2,3)
               ZROT(1,3)=AROT(3,1)
               ZROT(2,3)=AROT(3,2)
               ZROT(3,3)=AROT(3,3)
            ELSE
               DO 280 IHH=1,4
                  DO 270 JHH=1,4
                     ZROT(IHH,JHH)=AROT(IHH,JHH)
  270             CONTINUE
  280          CONTINUE
               CALL RINV(4,ZROT,4,RRR,IFAIL)
            ENDIF
            YROTMI=100000.
            YROTMA=-100000.
            DO 290 M=1,LMN
               IF(ISWIT(10).EQ.100)THEN
                  CXROT(M)=AROT(1,1)*CX(M)+ AROT(1,2)*CY(M)+AROT(1,3)*
     +            CZ(M)
                  CYROT(M)=AROT(2,1)*CX(M)+ AROT(2,2)*CY(M)+AROT(2,3)*
     +            CZ(M)
                  CZROT(M)=AROT(3,1)*CX(M)+ AROT(3,2)*CY(M)+AROT(3,3)*
     +            CZ(M)
               ELSE
                  CXROT(M)=AROT(1,4)+AROT(1,1)*CX(M)+ AROT(1,2)*CY(M)+
     +            AROT(1,3)*CZ(M)
                  CYROT(M)=AROT(2,4)+AROT(2,1)*CX(M)+ AROT(2,2)*CY(M)+
     +            AROT(2,3)*CZ(M)
                  CZROT(M)=AROT(3,4)+AROT(3,1)*CX(M)+ AROT(3,2)*CY(M)+
     +            AROT(3,3)*CZ(M)
               ENDIF
               IF(CYROT(M).LT.YROTMI)YROTMI=CYROT(M)
               IF(CYROT(M).GT.YROTMA)YROTMA=CYROT(M)
  290       CONTINUE
            IF(LINFIL.EQ.1)THEN
               RINULI=.01
            ELSEIF(LINFIL.EQ.2)THEN
               RINULI=.1
            ELSEIF(LINFIL.EQ.3)THEN
               RINULI=.05
            ELSEIF(LINFIL.EQ.4)THEN
               RINULI=.05
            ELSEIF(LINFIL.EQ.5)THEN
               RINULI=.02
            ELSEIF(LINFIL.EQ.6)THEN
               RINULI=.01
            ELSEIF(LINFIL.EQ.7)THEN
               RINULI=.005
            ENDIF
            NYROTM=(YROTMA-YROTMI)/RINULI
            YROTST=YROTMI
            DO 390 MM=1,NYROTM-1
               YROTST=YROTST+RINULI
               JFK=0
               DO 300 MMI=1,LMN-1,2
                  IF(CXROT(MMI).EQ.CXROT(MMI+1))THEN
                     IF(MIN(CYROT(MMI),CYROT(MMI+1)).LT. YROTST.AND.MAX
     +               (CYROT(MMI),CYROT(MMI+1)) .GT.YROTST)THEN
                        JFK=JFK+1
                        PXINT(JFK)=CXROT(MMI)
                        IF(CZROT(MMI).EQ.CZROT(MMI+1))THEN
                           PZINT(JFK)=CZROT(MMI)
                        ELSEIF(CYROT(MMI).NE.CYROT(MMI+1))THEN
                           AAAZ=(CYROT(MMI)-CYROT(MMI+1))/ (CZROT(MMI)-
     +                     CZROT(MMI+1))
                           BBBZ=(CZROT(MMI)*CYROT(MMI+1)- CYROT(MMI)*
     +                     CZROT(MMI+1))/ (CZROT(MMI)-CZROT(MMI+1))
                           PZINT(JFK)=(YROTST-BBBZ)/AAAZ
                        ENDIF
                     ENDIF
                  ELSEIF(CYROT(MMI).NE.CYROT(MMI+1))THEN
                     IF(MIN(CYROT(MMI),CYROT(MMI+1)).LT. YROTST.AND.MAX
     +               (CYROT(MMI),CYROT(MMI+1)) .GT.YROTST)THEN
                        AAA=(CYROT(MMI)-CYROT(MMI+1))/ (CXROT(MMI)-
     +                  CXROT(MMI+1))
                        BBB=(CXROT(MMI)*CYROT(MMI+1)- CYROT(MMI)*
     +                  CXROT(MMI+1))/ (CXROT(MMI)-CXROT(MMI+1))
                        JFK=JFK+1
                        PXINT(JFK)=(YROTST-BBB)/AAA
                        IF(CZROT(MMI).EQ.CZROT(MMI+1))THEN
                           PZINT(JFK)=CZROT(MMI)
                        ELSEIF(CYROT(MMI).NE.CYROT(MMI+1))THEN
                           AAAZ=(CYROT(MMI)-CYROT(MMI+1))/ (CZROT(MMI)-
     +                     CZROT(MMI+1))
                           BBBZ=(CZROT(MMI)*CYROT(MMI+1)- CYROT(MMI)*
     +                     CZROT(MMI+1))/ (CZROT(MMI)-CZROT(MMI+1))
                           PZINT(JFK)=(YROTST-BBBZ)/AAAZ
                        ENDIF
                     ENDIF
                  ENDIF
  300          CONTINUE
               DO 320 JM=1,JFK-1
                  DO 310 KM=JM+1,JFK
                     IF(PXINT(JM).LT.PXINT(KM))THEN
                        TMPM=PXINT(JM)
                        TMPZ=PZINT(JM)
                        PXINT(JM)=PXINT(KM)
                        PZINT(JM)=PZINT(KM)
                        PXINT(KM)=TMPM
                        PZINT(KM)=TMPZ
                     ENDIF
  310             CONTINUE
  320          CONTINUE
               RJFK=JFK*.5
               IRJFK=RJFK
*              IF((RJFK-IRJFK).GT..1)PRINT *,'Odd !'
               DO 330 MR=1,JFK
                  IF(ISWIT(10).EQ.100)THEN
                     XPINT(MR)=ZROT(1,1)*PXINT(MR)+ ZROT(1,2)*YROTST+
     +               ZROT(1,3)*CZROT(1)
                     YPINT(MR)=ZROT(2,1)*PXINT(MR)+ ZROT(2,2)*YROTST+
     +               ZROT(2,3)*CZROT(1)
                     ZPINT(MR)=ZROT(3,1)*PXINT(MR)+ ZROT(3,2)*YROTST+
     +               ZROT(3,3)*CZROT(1)
                  ELSE
                     XPINT(MR)=ZROT(1,4)+ZROT(1,1)*PXINT(MR)+ ZROT(1,2)
     +               *YROTST+ZROT(1,3)*PZINT(MR)
                     YPINT(MR)=ZROT(2,4)+ZROT(2,1)*PXINT(MR)+ ZROT(2,2)
     +               *YROTST+ZROT(2,3)*PZINT(MR)
                     ZPINT(MR)=ZROT(3,4)+ZROT(3,1)*PXINT(MR)+ ZROT(3,2)
     +               *YROTST+ZROT(3,3)*PZINT(MR)
                  ENDIF
  330          CONTINUE
               IF(LINFIL.GT.1)THEN
                  LINCOL=IBITS(LINATT,16,8)
                  CALL GDSHAD(LINCOL,APROSC)
               ENDIF
               DO 380 MZ=1,JFK,2
                  P1(1)=XPINT(MZ)
                  P1(2)=YPINT(MZ)
                  P1(3)=ZPINT(MZ)
                  P2(1)=XPINT(MZ+1)
                  P2(2)=YPINT(MZ+1)
                  P2(3)=ZPINT(MZ+1)
                  CALL CGVEDG(NT,P1,P2,IVIS)
                  IF (IVIS .LT. 0) GOTO 370
                  IF (NFACE .EQ. 0) GOTO 360
                  CALL CGVHED( Q(JCG+JXYZ1),IQ(JCG+JFA),IQ(JCG+JPFA),
     +            Q(JCG+ JDFA),IQ(JCG+JTRE),IQ(JCG+JSTA), Q(JCG+JALE),
     +            Q(JCG+JARI))
*           D R A W   E D G E
                  IF (NPART) 370 ,360 ,340
  340             DO 350  I=1,NPART+1
                     IF (I .EQ. 1) T1 = 0.
                     IF (I .NE. 1) T1 = TEND(I-1)
                     IF (I .NE. NPART+1) T2 = TSTRT(I)
                     IF (I .EQ. NPART+1) T2 = 1.
                     IF (T2-T1 .LT. TDEL) GOTO 350
                     AX(1) = XA + T1*XDELT
                     AY(1) = YA + T1*YDELT
                     AX(2) = XA + T2*XDELT
                     AY(2) = YA + T2*YDELT
                     IF(LINFIL.GT.1.AND.IDVIEW.NE.0)THEN
                        LLLINA=LINATT
                        CALL MVBITS(LINCOL,0,8,LINATT,16)
                        CALL GDRAWV(AX,AY,2)
                        LINATT=LLLINA
                     ELSE
                        CALL GDRAWV(AX,AY,2)
                     ENDIF
  350             CONTINUE
                  GOTO 370
*           L I N E   I S   F U L L Y   V I S I B L E
  360             CONTINUE
                  AX(1) = AA(1)
                  AY(1) = AA(2)
                  AX(2) = BB(1)
                  AY(2) = BB(2)
                  IF(LINFIL.GT.1.AND.IDVIEW.NE.0)THEN
                     LLLINA=LINATT
                     CALL MVBITS(LINCOL,0,8,LINATT,16)
                     CALL GDRAWV(AX,AY,2)
                     LINATT=LLLINA
                  ELSE
                     CALL GDRAWV(AX,AY,2)
                  ENDIF
  370             CONTINUE
  380          CONTINUE
  390       CONTINUE
         ENDIF
      ENDIF
***      IOLEP=LEP
***      LEP=11
      IF(LINFIL.NE.0)THEN
         LTY=IBITS(LINATT,10,3)
         CALL MVBITS(6,0,3,LINATT,10)
         CALL GDRAWV(TESTX,TESTY,0)
         CALL MVBITS(LTY,0,3,LINATT,10)
***      LEP=IOLEP
      ENDIF
*SG
  400 RETURN
*
10000    FORMAT('  First word of WIRE less or equal 0 ')
10100    FORMAT(' N. of words of Hidden Structure less or equal 18 ')
*
      END
#endif

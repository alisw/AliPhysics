*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGFARE(NT,FACE,IVIS,ISHAPE)
************************************************************************
*                                                                      *
*     Name: CGFARE                                                     *
*     Author: S. Giani                           Date:    20.05.91     *
*                                                Revised:              *
*                                                                      *
*     Function: HIDDEN FACE REMOVAL algoritm                           *
*               and transformation to screen coordinates               *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input:  NT - number of transformation to screen coordinates      *
*             FACE - face                                              *
*                                                                      *
*     Output: IVIS - visibility flag                                   *
*                    1 - if visible face                               *
*                   -1 - if unvisible                                  *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
#include "geant321/cgcfac.inc"
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/cgctra.inc"
#include "geant321/gcspee.inc"
#include "geant321/gcmutr.inc"
***SG
      DIMENSION ACCMI1(6),ACCMI2(6),ACCMI3(6),
     +          ACCMA1(6),ACCMA2(6),ACCMA3(6)
      DIMENSION SMI(3),SMA(3),POMI(3),POMA(3)
      SAVE ACCMI1,ACCMI2,ACCMI3,ACCMA1,ACCMA2,ACCMA3
      SAVE POMI,POMA,ACCXT1,ACCXT2,ACCNT1,ACCNT2
***SG
      REAL      FACE(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION  T(4,3),A,B,C,S
#endif
#if defined(CERNLIB_SINGLE)
      REAL      T(4,3)
#endif
*-
      IVIS   = -1
      DO 120 I=1,4
         DO 110 J=1,3
            T(I,J) = TSCRN(I,J,NT)
  110    CONTINUE
  120 CONTINUE
*
***SG
**         HIDDEN FACE REMOVAL
*     Computing face scope and skipping if it's 'covered': this
*     can allow a great increase in speed and a great reduction
*     in number of memory words used.
*
         J = LCGFAC
         NTIM=NTIM+1
         SRFMI1 = FACE(J+KCGX1)
         SRFMI2 = FACE(J+KCGY1)
         SRFMI3 = FACE(J+KCGZ1)
         SRFMA1 = FACE(J+KCGX1)
         SRFMA2 = FACE(J+KCGY1)
         SRFMA3 = FACE(J+KCGZ1)
         NEDGE = FACE(KCGNE)
         DO 333 NE=1,NEDGE
            SRFMI1 = MIN(SRFMI1,FACE(J+KCGX1),FACE(J+KCGX2))
            SRFMI2 = MIN(SRFMI2,FACE(J+KCGY1),FACE(J+KCGY2))
            SRFMI3 = MIN(SRFMI3,FACE(J+KCGZ1),FACE(J+KCGZ2))
            SRFMA1 = MAX(SRFMA1,FACE(J+KCGX1),FACE(J+KCGX2))
            SRFMA2 = MAX(SRFMA2,FACE(J+KCGY1),FACE(J+KCGY2))
            SRFMA3 = MAX(SRFMA3,FACE(J+KCGZ1),FACE(J+KCGZ2))
            J = J + LCGEDG
  333    CONTINUE
*        If volume set limits
         IF(IPORLI.EQ.1)THEN
*         If no clipping or shifting or exploding mode
          IF(NAIN.EQ.0.AND.KSHIFT.EQ.0.AND.GBOOM.EQ.0)THEN
*          If volume created by cgbox
           IF(ISHAPE.LT.5.OR.ISHAPE.EQ.10)THEN
*           Set 'faces scope' for sublim faces created by cgbox
            ACCMI1(NTIM) = SRFMI1
            ACCMI2(NTIM) = SRFMI2
            ACCMI3(NTIM) = SRFMI3
            ACCMA1(NTIM) = SRFMA1
            ACCMA2(NTIM) = SRFMA2
            ACCMA3(NTIM) = SRFMA3
*           Set 'volume scope' for sublim faces of revolution
            POMI(1)=S1
            POMI(2)=S2
            POMI(3)=S3
            POMA(1)=SS1
            POMA(2)=SS2
            POMA(3)=SS3
            ACCXT2=SRAGMX
            ACCXT1=SRAGMN
            ACCNT1=RAINT1
            ACCNT2=RAINT2
*          If volume of revolution
           ELSE
*           Set 'faces scope' for sublim faces created by cgbox
            ACCMI1(4)=S1
            ACCMI1(6)=SS1
            ACCMI2(3)=SS2
            ACCMI2(5)=S2
            ACCMI3(1)=S3
            ACCMI3(2)=SS3
            ACCMA1(4)=S1
            ACCMA1(6)=SS1
            ACCMA2(3)=SS2
            ACCMA2(5)=S2
            ACCMA3(1)=S3
            ACCMA3(2)=SS3
*           Set 'volume scope' and 'radial scope' for sublim faces of revolution
            POMI(1)=S1
            POMI(2)=S2
            POMI(3)=S3
            POMA(1)=SS1
            POMA(2)=SS2
            POMA(3)=SS3
            ACCXT2=SRAGMX
            ACCXT1=SRAGMN
            ACCNT1=RAINT1
            ACCNT2=RAINT2
           ENDIF
*         If clipping or shifting or exploding mode on
          ELSE
*          Set 'volume scope' for all sublim faces, and 'radial scope' as
*          well for sublim faces of revolution
           POMI(1)=S1
           POMI(2)=S2
           POMI(3)=S3
           POMA(1)=SS1
           POMA(2)=SS2
           POMA(3)=SS3
           ACCXT2=SRAGMX
           ACCXT1=SRAGMN
           ACCNT1=RAINT1
           ACCNT2=RAINT2
          ENDIF
*        If volume is to be compared with limits
         ELSEIF(ISUBLI.EQ.1)THEN
*         If no clipping or shifting or exploding mode
          IF(NAIN.EQ.0.AND.KSHIFT.EQ.0.AND.GBOOM.EQ.0)THEN
*          If volume created by cgbox
           IF(ISHAPE.LT.5.OR.ISHAPE.EQ.10)THEN
*           Comparison face by face
            IF(NTIM.EQ.1)THEN
               IF(SRFMI3.GT.ACCMI3(NTIM).AND.SRFMA3
     +         .GT.ACCMA3(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.2)THEN
               IF(SRFMI3.LT.ACCMI3(NTIM).AND.SRFMA3
     +         .LT.ACCMA3(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.3)THEN
               IF(SRFMI2.LT.ACCMI2(NTIM).AND.SRFMA2
     +         .LT.ACCMA2(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.5)THEN
               IF(SRFMI2.GT.ACCMI2(NTIM).AND.SRFMA2
     +         .GT.ACCMA2(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.4)THEN
               IF(SRFMI1.GT.ACCMI1(NTIM).AND.SRFMA1
     +         .GT.ACCMA1(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.6)THEN
               IF(SRFMI1.LT.ACCMI1(NTIM).AND.SRFMA1
     +         .LT.ACCMA1(NTIM))GOTO 999
            ENDIF
            GOTO 888
*          If volume of revolution
           ELSE
*           Comparison with mother scopes
            SMI(1)=SRFMI1
            SMI(2)=SRFMI2
            SMI(3)=SRFMI3
            SMA(1)=SRFMA1
            SMA(2)=SRFMA2
            SMA(3)=SRFMA3
            EXTRA1=RMAX1
            EXTRA2=RMAX2
            ENTRA1=RMIN1
            ENTRA2=RMIN2
*           If mother was created by cgbox or if it was of revolution
            ISP=0
            DO 127 I=1,3
             SPMI=SMI(I)-POMI(I)
             SPMA=SMA(I)-POMA(I)
             ASPMI=ABS(SPMI)
             ASPMA=ABS(SPMA)
             SMIA=SMI(I)-SMA(I)
             ASMIA=ABS(SMIA)
             IF(SPMI.GE.-0.001.AND.SPMA.LE.0.001)THEN
              ISP=ISP+1
              IF(ASPMI.LE.0.001.OR.ASPMA.LE.0.001)THEN
               IF(ASMIA.LE.0.001)GOTO 888
              ENDIF
             ENDIF
  127       CONTINUE
            IF(ISP.EQ.3)THEN
*            If mother was of revolution
             IF(ACCXT2.NE.0)THEN
              IF(ISCOP.EQ.1.AND.(ISHAPE.EQ.11.OR.ISHAPE.EQ.12
     +        .OR.ISHAPE.EQ.7.OR.ISHAPE.EQ.8))THEN
              EXXT1=EXTRA1-ACCXT1
              EXXT2=EXTRA2-ACCXT2
              ENNT1=ENTRA1-ACCNT1
              ENNT2=ENTRA2-ACCNT2
              IF(EXXT1.LT.-0.001.AND.EXXT2.LT.-0.001.AND.
     +        ENNT1.GT.0.001.AND.ENNT2.GT.0.001)THEN
               GOTO 999
              ELSEIF(EXXT1.LT.-0.001.AND.EXXT2.LT.-0.001.AND.
     +        ACCNT1.LT.0.001.AND.ACCNT2.LT.0.001)THEN
               GOTO 999
              ELSE
               GOTO 888
              ENDIF
              ELSE
              DO 701 ITER=1,IPORNT
                EXXT1=EXTRA1-PORMAR(ITER)
                EXXT2=EXTRA2-PORMAR(ITER)
                AEXXT1=ABS(EXXT1)
                AEXXT2=ABS(EXXT2)
                ENNT1=ENTRA1-PORMIR(ITER)
                ENNT2=ENTRA2-PORMIR(ITER)
                AENNT1=ABS(ENNT1)
                AENNT2=ABS(ENNT2)
                IF(AEXXT1.LT.0.001.OR.AEXXT2.LT.0.001)GOTO 888
                IF(AENNT1.LT.0.001.OR.AENNT2.LT.0.001)THEN
                  IF(PORMIR(ITER).NE.0.)GOTO 888
                ENDIF
 701          CONTINUE
              ENDIF
             ENDIF
             GOTO 999
            ELSE
             GOTO 888
            ENDIF
           ENDIF
*         If clipping or shifting or exploding mode on
          ELSE
*          Get scopes of the daughter (of each kind)
           SMI(1)=SRFMI1
           SMI(2)=SRFMI2
           SMI(3)=SRFMI3
           SMA(1)=SRFMA1
           SMA(2)=SRFMA2
           SMA(3)=SRFMA3
           EXTRA1=RMAX1
           EXTRA2=RMAX2
           ENTRA1=RMIN1
           ENTRA2=RMIN2
* If mother was clipped, check relative position of daughter and clipping
*   volumes; only if they don't interact, hidden face removal can work.
           DO 111 IJ=1,JPORJJ
            IFVFUN=0
            DO 301 J=1,3
             PMISMA=CLIPMI(J+3*IJ-3)-SMA(J)
             SMIPMA=SMI(J)-CLIPMA(J+3*IJ-3)
             APMISM=ABS(PMISMA)
             ASMIPM=ABS(SMIPMA)
             SMASMI=SMA(J)-SMI(J)
             ASMASM=ABS(SMASMI)
             IF(PMISMA.GE.-0.001.OR.
     +       SMIPMA.GE.-0.001)THEN
               IFVFUN=1
               IF(APMISM.LT.0.001.OR.
     +         ASMIPM.LT.0.001)THEN
                IF(ASMASM.LT.0.0001)GOTO 888
               ENDIF
               GO TO 302
             ENDIF
  301       CONTINUE
  302       CONTINUE
           IF(IFVFUN.EQ.0.AND.NAIN.NE.3)GO TO 888
  111      CONTINUE
*          If mother was created by cgbox or if it was of revolution
           ISP=0
           DO 128 I=1,3
            SPMI=SMI(I)-POMI(I)
            SPMA=SMA(I)-POMA(I)
            ASPMI=ABS(SPMI)
            ASPMA=ABS(SPMA)
            SMIA=SMI(I)-SMA(I)
            ASMIA=ABS(SMIA)
            IF(SPMI.GE.-0.001.AND.SPMA.LE.0.001)THEN
             ISP=ISP+1
             IF(ASPMI.LE.0.001.OR.ASPMA.LE.0.001)THEN
              IF(ASMIA.LE.0.001)GOTO 888
             ENDIF
            ENDIF
  128      CONTINUE
           IF(ISP.EQ.3)THEN
*           If mother was of revolution
            IF(ACCXT2.NE.0)THEN
              IF(ISCOP.EQ.1.AND.(ISHAPE.EQ.11.OR.ISHAPE.EQ.12
     +        .OR.ISHAPE.EQ.7.OR.ISHAPE.EQ.8))THEN
              EXXT1=EXTRA1-ACCXT1
              EXXT2=EXTRA2-ACCXT2
              ENNT1=ENTRA1-ACCNT1
              ENNT2=ENTRA2-ACCNT2
              IF(EXXT1.LT.-0.001.AND.EXXT2.LT.-0.001.AND.
     +        ENNT1.GT.0.001.AND.ENNT2.GT.0.001)THEN
               GOTO 999
              ELSEIF(EXXT1.LT.-0.001.AND.EXXT2.LT.-0.001.AND.
     +        ACCNT1.LT.0.001.AND.ACCNT2.LT.0.001)THEN
               GOTO 999
              ELSE
               GOTO 888
              ENDIF
              ELSE
              DO 702 ITER=1,IPORNT
                EXXT1=EXTRA1-PORMAR(ITER)
                EXXT2=EXTRA2-PORMAR(ITER)
                AEXXT1=ABS(EXXT1)
                AEXXT2=ABS(EXXT2)
                ENNT1=ENTRA1-PORMIR(ITER)
                ENNT2=ENTRA2-PORMIR(ITER)
                AENNT1=ABS(ENNT1)
                AENNT2=ABS(ENNT2)
                IF(AEXXT1.LT.0.001.OR.AEXXT2.LT.0.001)GOTO 888
                IF(AENNT1.LT.0.001.OR.AENNT2.LT.0.001)THEN
                  IF(PORMIR(ITER).NE.0.)GOTO 888
                ENDIF
 702          CONTINUE
              ENDIF
            ENDIF
            IF(ISCOP.EQ.1)THEN
             IF((ISHAPE.GT.1.AND.ISHAPE.LT.5).OR.ISHAPE.EQ.10)
     +       GOTO 888
            ENDIF
            GOTO 999
           ELSE
            GOTO 888
           ENDIF
          ENDIF
         ENDIF
 888  CONTINUE
*
***SG
*
      C      = (T(2,1)*T(3,2) - T(3,1)*T(2,2))*FACE(KCGAA) +
     +         (T(3,1)*T(1,2) - T(1,1)*T(3,2))*FACE(KCGBB) +
     +         (T(1,1)*T(2,2) - T(2,1)*T(1,2))*FACE(KCGCC)
      IF (C .LE. 0.)     GOTO 999
      B      = (T(2,3)*T(3,1) - T(3,3)*T(2,1))*FACE(KCGAA) +
     +         (T(3,3)*T(1,1) - T(1,3)*T(3,1))*FACE(KCGBB) +
     +         (T(1,3)*T(2,1) - T(2,3)*T(1,1))*FACE(KCGCC)
      A      = (T(2,2)*T(3,3) - T(3,2)*T(2,3))*FACE(KCGAA) +
     +         (T(3,2)*T(1,3) - T(1,2)*T(3,3))*FACE(KCGBB) +
     +         (T(1,2)*T(2,3) - T(2,2)*T(1,3))*FACE(KCGCC)
      S      = 1./SQRT(A*A+B*B+C*C)
      AABCD(1) = A*S
      AABCD(2) = B*S
      AABCD(3) = C*S
*
      F1(KCGAF) = FACE(KCGAF)
      F1(KCGNE) = FACE(KCGNE)
      F1(KCGAA) = 0.
      F1(KCGBB) = 0.
      F1(KCGCC) = 1.
      F1(KCGDD) = 0.
      F1(KCGNE) = FACE(KCGNE)
*
**           T R A S F E R   P O I N T   C O O R D I N A T E S
*
      NEDGE  = FACE(KCGNE)
      IF (LCGFAC+NEDGE*LCGEDG .GT. LABC)
     +       PRINT *, ' Problem in CGFVIS: no space'
      XGRAV  = 0.
      YGRAV  = 0.
      ZGRAV  = 0.
      J      = LCGFAC
      DO 200 NE=1,NEDGE
         F1(J+KCGAE) = FACE(J+KCGAE)
         X = FACE(J+KCGX1)
         Y = FACE(J+KCGY1)
         Z = FACE(J+KCGZ1)
         F1(J+KCGX1) = T(1,1)*X + T(2,1)*Y + T(3,1)*Z + T(4,1)
         F1(J+KCGY1) = T(1,2)*X + T(2,2)*Y + T(3,2)*Z + T(4,2)
         F1(J+KCGZ1) = T(1,3)*X + T(2,3)*Y + T(3,3)*Z + T(4,3)
         X = FACE(J+KCGX2)
         Y = FACE(J+KCGY2)
         Z = FACE(J+KCGZ2)
         F1(J+KCGX2) = T(1,1)*X + T(2,1)*Y + T(3,1)*Z + T(4,1)
         F1(J+KCGY2) = T(1,2)*X + T(2,2)*Y + T(3,2)*Z + T(4,2)
         F1(J+KCGZ2) = T(1,3)*X + T(2,3)*Y + T(3,3)*Z + T(4,3)
         XGRAV = XGRAV + F1(J+KCGX1) + F1(J+KCGX2)
         YGRAV = YGRAV + F1(J+KCGY1) + F1(J+KCGY2)
         ZGRAV = ZGRAV + F1(J+KCGZ1) + F1(J+KCGZ2)
         J = J + LCGEDG
  200 CONTINUE
      XFACT = 1./(2.*NEDGE)
      XGRAV = XGRAV * XFACT
      YGRAV = YGRAV * XFACT
      ZGRAV = ZGRAV * XFACT
      AABCD(4) =-(AABCD(1)*XGRAV + AABCD(2)*YGRAV + AABCD(3)*ZGRAV)
*
**           F I N D   F A C E   M I N - M A X
*
      J      = LCGFAC
      RFMIN(1) = F1(J+KCGX1)
      RFMIN(2) = F1(J+KCGY1)
      RFMIN(3) = F1(J+KCGZ1)
      RFMAX(1) = F1(J+KCGX1)
      RFMAX(2) = F1(J+KCGY1)
      RFMAX(3) = F1(J+KCGZ1)
      DO 300 NE=1,NEDGE
         RFMIN(1) = MIN(RFMIN(1),F1(J+KCGX1),F1(J+KCGX2))
         RFMIN(2) = MIN(RFMIN(2),F1(J+KCGY1),F1(J+KCGY2))
         RFMIN(3) = MIN(RFMIN(3),F1(J+KCGZ1),F1(J+KCGZ2))
         RFMAX(1) = MAX(RFMAX(1),F1(J+KCGX1),F1(J+KCGX2))
         RFMAX(2) = MAX(RFMAX(2),F1(J+KCGY1),F1(J+KCGY2))
         RFMAX(3) = MAX(RFMAX(3),F1(J+KCGZ1),F1(J+KCGZ2))
         F1(J+KCGZ1) = 0.
         F1(J+KCGZ2) = 0.
         J = J + LCGEDG
  300 CONTINUE
      DRFACE(1) =-RFMAX(1)
      DRFACE(2) =-RFMAX(2)
      DRFACE(3) = RFMIN(1)
      DRFACE(4) = RFMIN(2)
      DRFACE(5) = RFMIN(3)
*
**          C O M P U T E    F A C E    V I S I B L E    A R E A
*
      J      = LCGFAC
      S      = 0.
      DLMAX  = 0.
      DO 400 NE=1,NEDGE
         S = S + F1(J+KCGX1)*F1(J+KCGY2) - F1(J+KCGX2)*F1(J+KCGY1)
         DL = ABS(F1(J+KCGX2)-F1(J+KCGX1)) + ABS(F1(J+KCGY2)-F1(J+
     +   KCGY1))
         IF (DLMAX .LT. DL) DLMAX = DL
         J = J + LCGEDG
  400 CONTINUE
      IF (DLMAX .LT. EESCR)     GOTO 999
      IF (S .GT. DLMAX*EESCR)   IVIS = 1
*
  999 RETURN
      END

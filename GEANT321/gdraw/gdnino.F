*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 28/03/94  01.30.59  by  S.Giani
*-- Author :
      SUBROUTINE GDNINO
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       RAY-TRACING                                              *
C.    *                                                                *
C.    *       Tracking of light rays.                                  *
C.    *                                                                *
C.    *       Author: S.Giani.                                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#if defined(CERNLIB_DEBUG)
#include "geant321/gcunit.inc"
#endif
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
********************************************************************************
#include "geant321/gcfdim.inc"
#include "geant321/gcrayt.inc"
#include "geant321/gcpixe.inc"
      DIMENSION CCVV(3),DDVV(2)
********************************************************************************
C.
C.    ------------------------------------------------------------------
C.
      STEP = BIG
      IF((BOFLAG.NE.1.).and.(VECT(1).le.XCUT.or.VECT(2).le.YCUT
     +.or.VECT(3).le.ZCUT))THEN
       CALL GTNEXT
      ELSE
       IGNEXT=1
       SNEXT=0.0001
       nlevin=nlevel
      ENDIF
      IF (IGNEXT.NE.0) THEN
         STEP = SNEXT + PREC
         step=step+.004
***         if(iombra.eq.9)step=step+.001
         INWVOL= 2
         LMEC(1)=1
*
* *** Linear transport
*
         DO 10 I = 1,3
            VECTMP  = VECT(I) +STEP*VECT(I+3)
            IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
               IF(VECT(I+3).NE.0.) THEN
                  VECTMP = VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*
     +            EPSMAC
               ENDIF
            ENDIF
********************************************************************************
            IF(BOFLAG.NE.1)THEN
              VECT(I) = VECTMP
            ELSE
              INWVOL=0
            ENDIF
   10    CONTINUE
       NOFLAG=0
       IF(BOFLAG.NE.1)THEN
       IF(VECT(1).GT.XCUT.AND.VECT(2).GT.YCUT.AND.
     +    VECT(3).GT.ZCUT)THEN
          NOFLAG=1
          CCXX(1)=1.
          CCXX(2)=1.
          CCXX(3)=0.
          CCXX(4)=0.
          CCYY(1)=1.
          CCYY(2)=0.
          CCYY(3)=1.
          CCYY(4)=0.
          CCZZ(1)=1.
          CCZZ(2)=0.
          CCZZ(3)=0.
          CCZZ(4)=1.
          DO 11 I=1,3
            CCVV(I)=VECT(I+3)
 11       CONTINUE
          STEP1=100000.
          DDXX=GVSAFE(VECT,-XCUT,CCXX,4)
          CALL GVLONG(CCVV,CCXX,4,DDVV,NTAU)
          IF(NTAU.GT.0)THEN
           IF(DDVV(1).GE.0.)THEN
             STEP1=DDVV(1)
           ENDIF
           APFLAG=(CCXX(2)*XCOSXS)+(CCXX(3)*YCOSYS)+
     +     (CCXX(4)*ZCOSZS)
           IXYFLA=1
          ENDIF
          DDYY=GVSAFE(VECT,-YCUT,CCYY,4)
          CALL GVLONG(CCVV,CCYY,4,DDVV,NTAU)
          IF(NTAU.GT.0)THEN
           IF(DDVV(1).LT.STEP1.AND.DDVV(1).GE.0.)THEN
             IXYFLA=2
             STEP1=DDVV(1)
             APFLAG=(CCYY(2)*XCOSXS)+(CCYY(3)*YCOSYS)+
     +       (CCYY(4)*ZCOSZS)
           ENDIF
          ENDIF
          DDZZ=GVSAFE(VECT,-ZCUT,CCZZ,4)
          CALL GVLONG(CCVV,CCZZ,4,DDVV,NTAU)
          IF(NTAU.GT.0)THEN
           IF(DDVV(1).LT.STEP1.AND.DDVV(1).GE.0.)THEN
             IXYFLA=3
             STEP1=DDVV(1)
             APFLAG=(CCZZ(2)*XCOSXS)+(CCZZ(3)*YCOSYS)+
     +       (CCZZ(4)*ZCOSZS)
           ENDIF
          ENDIF
          STEP1=STEP1+.004
          DO 12 I=1,3
            IF(BOFLAG.NE.1)VECT(I)=VECT(I)+STEP1*VECT(I+3)
 12       CONTINUE
          IF(IOMBRA.EQ.9)SLENG=SLENG+STEP1
          IF(IOMBRA.EQ.10)THEN
            XYZMOD=1./SQRT((XLPOS-VECT(1))**2+(YLPOS-VECT(2))**2+
     +                     (ZLPOS-VECT(3))**2)
            XLDIR=(XLPOS-VECT(1))*XYZMOD
            YLDIR=(YLPOS-VECT(2))*XYZMOD
            ZLDIR=(ZLPOS-VECT(3))*XYZMOD
            IF(IXYFLA.EQ.1)THEN
             APFLAH=(CCXX(2)*XLDIR)+(CCXX(3)*YLDIR)+(CCXX(4)*ZLDIR)
*             if(apflah.lt.0.)print *,APFLAH,'=apflah',IXYFLA,'=flag',
*     +       CCXX(2),XLDIR,CCXX(3),YLDIR,CCXX(4),ZLDIR
            ELSEIF(IXYFLA.EQ.2)THEN
             APFLAH=(CCYY(2)*XLDIR)+(CCYY(3)*YLDIR)+(CCYY(4)*ZLDIR)
*             if(apflah.lt.0.)print *,APFLAH,'=apflah',IXYFLA,'=flag',
*     +       CCYY(2),XLDIR,CCYY(3),YLDIR,CCYY(4),ZLDIR
            ELSEIF(IXYFLA.EQ.3)THEN
             APFLAH=(CCZZ(2)*XLDIR)+(CCZZ(3)*YLDIR)+(CCZZ(4)*ZLDIR)
            ENDIF
          ENDIF
       ENDIF
       ELSE
        IF(IOMBRA.EQ.10)THEN
          IOMBRA=9
          XYZMOD=1./SQRT((XLPOS-VECT(1))**2+(YLPOS-VECT(2))**2+
     +                   (ZLPOS-VECT(3))**2)
          VECT(4)=(XLPOS-VECT(1))*XYZMOD
          VECT(5)=(YLPOS-VECT(2))*XYZMOD
          VECT(6)=(ZLPOS-VECT(3))*XYZMOD
          XLDIR=VECT(4)
          YLDIR=VECT(5)
          ZLDIR=VECT(6)
          SLENG=0.
        ELSE
          NSTEP=10000
          ISTOP=0
          IF(IOMBRA.EQ.9)IOMBRA=10
        ENDIF
       ENDIF
      ELSE
          NSTEP=10000
          ISTOP=0
          INWVOL= 2
          LMEC(1)=1
          STEP=.1
          IF(IOMBRA.EQ.9)THEN
             IOMBRA=10
             PRINT *,'WARNING'
          ENDIF
         DO 20 I = 1,3
           IF(BOFLAG.NE.1)VECT(I)  = VECT(I) +STEP*VECT(I+3)
********************************************************************************
   20    CONTINUE
      ENDIF
      SLENG = SLENG +STEP
      END

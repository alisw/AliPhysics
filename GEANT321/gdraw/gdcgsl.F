*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCGSL(IVOLNA,ISHAPE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine allows computes the coefficients of the     *
C.    *       cut plane and the limits and array of the clipping       *
C.    *       volumes (boxes, cones, tubes, spheres).                  *
C.    *                                                                *
C.    *    ==>Called by : GDRAW                                        *
C.    *                                                                *
C.    *       Authors :  J.Salt ; S.Giani                              *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcgobj.inc"
*
*
*****SG
*
#include "geant321/gcspee.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcmutr.inc"
#include "geant321/pawc.inc"
*
      DIMENSION TMIN(3),TMAX(3),XZ(2,4)
      CHARACTER*4 NACA
*
      CALL UCTOH('PERS',IPERS,4,4)
      IF(IHOLE.NE.0)THEN
*
*   Clipping Volumes Creation
*
*   Look for volume to be clipped
*
         NAIN=0
         JJJ=0
         ICUBE=0
         CALL UHTOC(IVOLNA,4,NACA,4)
         DO 30  III=1,MULTRA
            INV=0
            IWILDC=INDEX(GNNVV(III),'*')
            IF(IWILDC.EQ.0) THEN
               CALL UCTOH(GNNVV(III),INV,4,4)
            ELSEIF(IWILDC.EQ.1.AND.JJJ.LT.2) THEN
               INV=IVOLNA
            ELSEIF(GNNVV(III)(1:IWILDC-1).EQ.
     +                NACA(1:IWILDC-1)) THEN
               INV=IVOLNA
            ENDIF
            IF(INV.EQ.IVOLNA)THEN
*
*   If you find it, compute number of times it's to be clipped
*   and set parameters of relative clipping shapes
*
               JJJ=JJJ+1
               IF(JJJ.EQ.3)THEN
                  WRITE(CHMAIL,10100)
                  CALL GMAIL(0,0)
                  GOTO 60
               ENDIF
               IF(GNASH(III).EQ.'BOX ')THEN
                  DBX=GXMAX(III)-GXMIN(III)
                  DBY=GYMAX(III)-GYMIN(III)
                  DBZ=GZMAX(III)-GZMIN(III)
*
                  IF(JJJ.EQ.2)THEN
                     ICUBE=ICUBE+1
                     CALL CGBRIK(DBX,DBY,DBZ,300,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     CALL CGSHIF(GXMIN(III),GYMIN(III),GZMIN(III),
     +               Q(ICLIP2))
                     CALL CGMNMX(Q(ICLIP2),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP2))
                     ENDIF
                  ELSE
                     ICUBE=ICUBE+1
                     CALL CGBRIK(DBX,DBY,DBZ,300,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     CALL CGSHIF(GXMIN(III),GYMIN(III),GZMIN(III),
     +               Q(ICLIP1))
                     CALL CGMNMX(Q(ICLIP1),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP1))
                     ENDIF
                  ENDIF
               ELSE IF (GNASH(III).EQ.'TUBE') THEN
                  RMIN1=0
                  RMAX1=GXMIN(III)
                  Z2=GXMAX(III)
                  RMIN2=RMIN1
                  RMAX2=RMAX1
                  XZ(1,1)=RMIN1
                  XZ(2,1)=-Z2
                  XZ(1,2)=RMAX1
                  XZ(2,2)=-Z2
                  XZ(1,3)=RMAX2
                  XZ(2,3)=Z2
                  XZ(1,4)=RMIN2
                  XZ(2,4)=Z2
                  ANG1=0.
                  ANG2=360.
                  NANG=30
*
                  IF(JJJ.EQ.2)THEN
                     CALL CGZREV(XZ,ANG1,ANG2,NANG,16000,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     S1=GYMIN(III)
                     S2=GYMAX(III)
                     S3=GZMIN(III)
                     CALL CGSHIF(S1,S2,S3,Q(ICLIP2))
                     CALL CGMNMX(Q(ICLIP2),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP2))
                     ENDIF
                  ELSE
                     CALL CGZREV(XZ,ANG1,ANG2,NANG,16000,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     S1=GYMIN(III)
                     S2=GYMAX(III)
                     S3=GZMIN(III)
                     CALL CGSHIF(S1,S2,S3,Q(ICLIP1))
                     CALL CGMNMX(Q(ICLIP1),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP1))
                     ENDIF
                  ENDIF
               ELSE IF (GNASH(III).EQ.'SPHE') THEN
                  R=GXMIN(III)
                  NLAT=15
                  NLON=15
*
                  IF(JJJ.EQ.2)THEN
                     CALL CGSPHE(R,NLAT,NLON,16000,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     S1=GXMAX(III)
                     S2=GYMIN(III)
                     S3=GYMAX(III)
                     CALL CGSHIF(S1,S2,S3,Q(ICLIP2))
                     CALL CGMNMX(Q(ICLIP2),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP2))
                     ENDIF
                  ELSE
                     CALL CGSPHE(R,NLAT,NLON,16000,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     S1=GXMAX(III)
                     S2=GYMIN(III)
                     S3=GYMAX(III)
                     CALL CGSHIF(S1,S2,S3,Q(ICLIP1))
                     CALL CGMNMX(Q(ICLIP1),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP1))
                     ENDIF
                  ENDIF
               ELSE IF (GNASH(III).EQ.'CONE') THEN
                  RMIN1=0.
                  RMAX1=GXMIN(III)
                  RMIN2=0.
                  RMAX2=GXMAX(III)
                  Z2=GYMIN(III)
                  XZ(1,1)=RMIN1
                  XZ(2,1)=-Z2
                  XZ(1,2)=RMAX1
                  XZ(2,2)=-Z2
                  XZ(1,3)=RMAX2
                  XZ(2,3)=Z2
                  XZ(1,4)=RMIN2
                  XZ(2,4)=Z2
                  ANG1=0.
                  ANG2=360.
                  NANG=30
                  IF(JJJ.EQ.2)THEN
                     CALL CGZREV(XZ,ANG1,ANG2,NANG,16000,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     CALL CGCEV(-1,Q(ICLIP2))
                     S1=GYMAX(III)
                     S2=GZMIN(III)
                     S3=GZMAX(III)
                     CALL CGSHIF(S1,S2,S3,Q(ICLIP2))
                     CALL CGMNMX(Q(ICLIP2),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP2))
                     ENDIF
                  ELSE
                     CALL CGZREV(XZ,ANG1,ANG2,NANG,16000,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     CALL CGCEV(-1,Q(ICLIP1))
                     S1=GYMAX(III)
                     S2=GZMIN(III)
                     S3=GZMAX(III)
                     CALL CGSHIF(S1,S2,S3,Q(ICLIP1))
                     CALL CGMNMX(Q(ICLIP1),TMIN,TMAX)
*   Perspective view
                     IF (IPRJ.EQ.IPERS) THEN
                        CALL CGPERS(Q(ICLIP1))
                     ENDIF
                  ENDIF
               ENDIF
*          IF(CGERR.LE.0)THEN
*             CALL GDCGER(CGERR)
*             IF(KCGST.EQ.-2)GO TO 999
*             IF(KCGST.EQ.-3)THEN
*                WRITE(CHMAIL,10100)
*                CALL GMAIL(0,0)
*                GO TO 999
*             ENDIF
*          ENDIF
*
*    Compute scope for each clipping volume
*
               DO 10 K=1,3
                  KKK=K+3*JJJ-3
                  BMIN(KKK)=TMIN(K)
                  BMAX(KKK)=TMAX(K)
   10          CONTINUE
               IF(IPORLI.EQ.1)THEN
                  DO 20 KJ=1,3
                     KKKJ=KJ+3*JJJ-3
                     CLIPMI(KKKJ)=TMIN(KJ)
                     CLIPMA(KKKJ)=TMAX(KJ)
   20             CONTINUE
               ENDIF
            ENDIF
   30    CONTINUE
*
*    If volume is not to be clipped
*
         IF(JJJ.EQ.0)THEN
            IF(IPORLI.EQ.1)THEN
               DO 40 KJ=1,6
                  CLIPMI(KJ)=-10000
                  CLIPMA(KJ)=-9999
   40          CONTINUE
               JPORJJ=1
            ENDIF
            NAIN=1
            GOTO 60
         ELSE
            IF(IPORLI.EQ.1)THEN
               JPORJJ=JJJ
            ENDIF
            NAIN=2
            ISA=0
            DO 50  J=1,6
               IF(BMIN(J).EQ.CLIPMI(J).AND.BMAX(J).EQ.CLIPMA(J))THEN
                  ISA=ISA+1
               ENDIF
   50       CONTINUE
            IF(ISA.EQ.6)NAIN=3
         ENDIF
*
*
*****SG
*
*
      ELSE
*
*   Slicing with a plane
*
         IF(ICUT.EQ.0) GO TO 999
         IF(ICUT.EQ.1)THEN
            CTHETA=90.
            CPHI=180.
         ELSE IF(ICUT.EQ.2)THEN
            CTHETA=90.
            CPHI=270.
         ELSE IF(ICUT.EQ.3)THEN
            CTHETA=180.
            CPHI=0.
         ELSE
            WRITE(CHMAIL,10000)
            CALL GMAIL(0,0)
            GO TO 999
         ENDIF
         ATH=DEGRAD*CTHETA
         APH=DEGRAD*CPHI
         ABCD(1)=SIN(ATH)*COS(APH)
         ABCD(2)=SIN(ATH)*SIN(APH)
         ABCD(3)=COS(ATH)
         ABCD(4)=DCUT
      ENDIF
   60 CONTINUE
*
10000 FORMAT(' CUT Index not implemented')
*10100 FORMAT(' Check Clipping Box Parameters ')
10100 FORMAT(' Please, reset CVOL mode. ')
*
  999 END

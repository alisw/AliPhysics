*
* $Id$
*
* $Log$
* Revision 1.2  1996/02/27 10:02:05  ravndal
* Drawing of PCON's optimized for 'HIDE ON'
*
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 13/12/94  17.13.38  by  S.Giani
*-- Author :
*
      SUBROUTINE GDCGOB(ITASK,ISHAPE,PAR,NOBJ,NWWS,IVOLNA,
     +LSTEP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Make the CG-Object with shape ISHAPE of parameters PAR   *
C.    *       with the same logic as GDRAWS. 1992                      *
C.    *                                                                *
C.    *     Input Parameters :                                         *
C.    *                                                                *
C.    *     ITASK:   Number for indicating the task to be performed    *
C.    *                                                                *
C.    *                                                                *
C.    *           = 0      Counting task                               *
C.    *           = 1      Slicing + Counting                          *
C.    *           = 2      Clipping + Counting                         *
C.    *           = 3      Insert into the H.S. + Convert to Wire      *
C.    *           = 4      Slicing + Insert into the H.S. + Convert    *
C.    *                    to Wire                                     *
C.    *           = 5      Clipping + Insert into the H.S. + Convert   *
C.    *                    to Wire                                     *
C.    *                                                                *
C.    * SHAPE     SHAPE    SHAPE                                       *
C.    * NUMBER    TYPE     PARAMETERS                                  *
C.    * -------------------------------------------------------------- *
C.    *                                                                *
C.    *   1       BOX      DX,DY,DZ                                    *
C.    *   2       TRD1     DX1,DX2,DY,DZ                               *
C.    *   3       TRD2     DX1,DX2,DY1,DY2,DZ                          *
C.    *   4       TRAP     DZ,TX,TY,H1,BL1,TL1,TTH1,H2,BL2,TL2,TTH2    *
C.    *                                                                *
C.    *   5       TUBE     RMIN,RMAX,DZ                                *
C.    *   6       TUBS     RMIN,RMAX,DZ,PHIMIN,PHIMAX                  *
C.    *   7       CONE     DZ,RMIN1,RMAX1,RMIN2,RMAX2                  *
C.    *   8       CONS     DZ,RMIN1,RMAX1,RMIN2,RMAX2,PHIMIN,PHIMAX    *
C.    *                                                                *
C.    *   9       SPHE     RMIN,RMAX,THEMIN,THEMAX,PHIMIN,PHIMAX       *
C.    *                                                                *
C.    *  10       PARA     DX,DY,DZ,TXY,TXZ,TYZ                        *
C.    *  11       PGON     PHIMIN,DPHI,NDIV,NZ,Z(1),RMIN(1),RMAX(1),...*
C.    *  12       PCON     PHIMIN,DPHI,NZ,Z(1),RMIN(1),RMAX(1),Z(2),...*
C.    *                                                                *
C.    * NOBJ = Counter of cg objects                                   *
C.    * NWWS = Size of Wire structure                                  *
C.    * IVOLNA = Name of volume                                        *
C.    * LSTEP = Number of CG objects forming each volume                *
C.    *                                                                *
C.    *    ==>Called by : GDRAW                                        *
C.    *       Author : P.Zanarini, J.Salt, S.Giani *********           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gcmutr.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
#include "geant321/gconsp.inc"
      SAVE NWPROD
      COMMON /QUEST/IQUEST(100)
      DIMENSION PAR(50),P(3,8)
*SG
      DIMENSION RRMIN(3),RRMAX(3)
      DIMENSION SLI1(4),SLI2(4),SPI1(4),SPI2(4)
*SG
      DIMENSION XZ(2,4),ZR(18),RMIR(18),RMAR(18),AMIRMA(18),AMARMA(18)
      DIMENSION T(4,3)
C.    ------------------------------------------------------------------
C.
**SG
      CALL UCTOH('PERS',IPERS,4,4)
      T(4,1)=0.
      T(4,2)=0.
      T(4,3)=0.
      LINSTY=IBITS(LINATT,10,3)
      IF(LINSTY.EQ.7)THEN
         APPROS=30.
      ELSE
         APPROS=15.
      ENDIF
      IF(ISUBLI.LT.IOLDSU)THEN
         PORGX=0
         PORGY=0
         PORGZ=0
         DO 10  J=1,15
            POX(J)=0
            POY(J)=0
            POZ(J)=0
   10    CONTINUE
      ENDIF
      IOLDSU=ISUBLI
*
*       LHC flag 'ON' (default)
*
*      CALL UCTOH('ON  ',LHIF,4,4)
*      IF(LEP.EQ.LHIF)THEN
*         VITE=1
*      ELSE
*         VITE=0
*      ENDIF
*
* Flag for GDCGHI resetted for each CG object
      ISG=0
      ICGP=0
      LINFIL=IBITS(LINATT,13,3)
**SG
      IVCLOS=0
      IVFUN=1
      IWORK=ITASK
      CALL UCTOH('ON  ',IFLH,4,4)
*JS
      IF (ISHAPE.EQ.1) THEN
C
C             BOX
C
         DX1=PAR(1)
         DY1=PAR(2)
         DX2=DX1
         DY2=DY1
         DZ=PAR(3)
         GO TO 20
C
      ELSEIF (ISHAPE.EQ.2) THEN
C
C             TRD1
C
         DX1=PAR(1)
         DX2=PAR(2)
         DY1=PAR(3)
         DY2=DY1
         DZ=PAR(4)
         GO TO 20
C
      ELSEIF (ISHAPE.EQ.3) THEN
C
C             TRD2
C
         DX1=PAR(1)
         DX2=PAR(2)
         DY1=PAR(3)
         DY2=PAR(4)
         DZ=PAR(5)
         GO TO 20
C
      ELSEIF (ISHAPE.EQ.4) THEN
C
C             TRAP
C
         DZ=PAR(1)
         TX=PAR(2)
         TY=PAR(3)
         H1=PAR(4)
         BL1=PAR(5)
         TL1=PAR(6)
         TTH1=PAR(7)
         H2=PAR(8)
         BL2=PAR(9)
         TL2=PAR(10)
         TTH2=PAR(11)
         GO TO 30
C
      ELSEIF (ISHAPE.EQ.5) THEN
C
C             TUBE
C
         AFINV=1./COS(PI/APPROS)
         FINV=ABS(AFINV)
         RMIN1=PAR(1)*FINV
         RMAX1=PAR(2)*FINV
         RMIN2=RMIN1
         RMAX2=RMAX1
         Z2=PAR(3)
*         Z1=-Z2
         PHIMIN=0.
         PHIMAX=360.
         IF((LINFIL.EQ.2.OR.LINFIL.EQ.3)
     +   .AND.RMIN1.NE.0)PHIMIN=5.
*SG
         ANG1=PHIMIN
         ANG2=PHIMAX
         AANG=ABS(ANG2-ANG1)
         AZLAT=AANG*APPROS
         ZLAT=AZLAT/360
         NANG=ZLAT
         IF(NANG.EQ.0)NANG=1
         AZ=ZLAT-NANG
         IF(AZ.GT..5)NANG=NANG+1
*SG
         GO TO 70
C
      ELSEIF (ISHAPE.EQ.6.OR.ISHAPE.EQ.29) THEN
C
C             TUBS
C
         AFINV=1./COS(PI/APPROS)
         FINV=ABS(AFINV)
         RMIN1=PAR(1)*FINV
         RMAX1=PAR(2)*FINV
         RMIN2=RMIN1
         RMAX2=RMAX1
         AZ2=PAR(3)
         Z2=ABS(AZ2)
*         Z1=-Z2
         PHIMIN=PAR(4)
         PHIMAX=PAR(5)
**SG
         ANG1=PHIMIN
         ANG2=PHIMAX
         AANG=ABS(ANG2-ANG1)
         AZLAT=AANG*APPROS
         ZLAT=AZLAT/360
         NANG=ZLAT
         IF(NANG.EQ.0)NANG=1
         AZ=ZLAT-NANG
         IF(AZ.GT..5)NANG=NANG+1
         IF(ISHAPE.EQ.29)NANG=APPROS
**SG
         GO TO 70
C
      ELSEIF (ISHAPE.EQ.7) THEN
C
C             CONE
C
         AFINV=1./COS(PI/APPROS)
         FINV=ABS(AFINV)
         RMIN1=PAR(2)*FINV
         RMAX1=PAR(3)*FINV
         RMIN2=PAR(4)*FINV
         RMAX2=PAR(5)*FINV
         Z2=PAR(1)
*         Z1=-Z2
         PHIMIN=0.
         PHIMAX=360.
         IF((LINFIL.EQ.2.OR.LINFIL.EQ.3).AND.RMIN2.NE.0)PHIMIN=5.
*SG
         ANG1=PHIMIN
         ANG2=PHIMAX
         AANG=ABS(ANG2-ANG1)
         AZLAT=AANG*APPROS
         ZLAT=AZLAT/360
         NANG=ZLAT
         IF(NANG.EQ.0)NANG=1
         AZ=ZLAT-NANG
         IF(AZ.GT..5)NANG=NANG+1
*SG
         GO TO 70
C
      ELSEIF (ISHAPE.EQ.8) THEN
C
C             CONS
C
         AFINV=1./COS(PI/APPROS)
         FINV=ABS(AFINV)
         RMIN1=PAR(2)*FINV
         RMAX1=PAR(3)*FINV
         RMIN2=PAR(4)*FINV
         RMAX2=PAR(5)*FINV
         Z2=PAR(1)
*         Z1=-Z2
         PHIMIN=PAR(6)
         PHIMAX=PAR(7)
**SG
         ANG1=PHIMIN
         ANG2=PHIMAX+.1
         AANG=ABS(ANG2-ANG1)
         AZLAT=AANG*APPROS
         ZLAT=AZLAT/360
         NANG=ZLAT
         IF(NANG.EQ.0)NANG=1
         AZ=ZLAT-NANG
         IF(AZ.GT..5)NANG=NANG+1
**SG
         GO TO 70
C
      ELSEIF (ISHAPE.EQ.9) THEN
C
C             SPHE
C
*         RMIN=PAR(1)
         RMAX=PAR(2)
         GO TO 120
C
      ELSEIF (ISHAPE.EQ.10) THEN
C
C             PARA
C
         DX=PAR(1)
         DY=PAR(2)
         DZ=PAR(3)
         TXY=PAR(4)
         TXZ=PAR(5)
         TYZ=PAR(6)
C
         TX=TXZ
         TY=TYZ
         H1=DY
         BL1=DX
         TL1=DX
         TTH1=TXY
         H2=DY
         BL2=DX
         TL2=DX
         TTH2=TXY
         GO TO 30
C
      ELSEIF (ISHAPE.EQ.11) THEN
C
C             PGON
C
         PHIMIN=PAR(1)
         IF((LINFIL.EQ.2.OR.LINFIL.EQ.3).AND.PAR(2).GT.359.)PAR(2)=359.
         PHIMAX=PHIMIN+PAR(2)
         NDIVAN=PAR(3)
         NZ=PAR(4)
C
C             Z(1)=PAR(5) - RMIN(1)=PAR(6) - RMAX(1)=PAR(7) - Z(2)=PAR(8
C
         GO TO 150
C
      ELSEIF (ISHAPE.EQ.12) THEN
C
C             PCON
C
         PHIMIN=PAR(1)
         IF((LINFIL.EQ.2.OR.LINFIL.EQ.3).AND.PAR(2).GT.359.)PAR(2)=359.
         PHIMAX=PHIMIN+PAR(2)
         NZ=PAR(3)
C
C             Z(1)=PAR(4) - RMIN(1)=PAR(5) - RMAX(1)=PAR(6) - Z(2)=PAR(7
C
         GO TO 230
      ELSE
         GO TO 999
      ENDIF
C
*      GO TO 150
C
   20 CONTINUE
C
C             Rectilinear shapes: BOX,TRD1,TRD2
C
      X1=0.
      Y1=0.
      X2=0.
      Y2=0.
      IF(DZ.LT.0.001)DZ=0.001
      Z1=-DZ
      Z2=DZ
C
C             Calculate the 8 vertex for rectilinear shapes
C
      IF(DX1.EQ.0.)DX1=0.0001
      IF(DY1.EQ.0.)DY1=0.0001
      IF(DX2.EQ.0.)DX2=0.0001
      IF(DY2.EQ.0.)DY2=0.0001
      P(1,1)=X1+DX1
      P(2,1)=Y1+DY1
      P(3,1)=Z1
      P(1,2)=X1-DX1
      P(2,2)=Y1+DY1
      P(3,2)=Z1
      P(1,3)=X1-DX1
      P(2,3)=Y1-DY1
      P(3,3)=Z1
      P(1,4)=X1+DX1
      P(2,4)=Y1-DY1
      P(3,4)=Z1
      P(1,5)=X2+DX2
      P(2,5)=Y2+DY2
      P(3,5)=Z2
      P(1,6)=X2-DX2
      P(2,6)=Y2+DY2
      P(3,6)=Z2
      P(1,7)=X2-DX2
      P(2,7)=Y2-DY2
      P(3,7)=Z2
      P(1,8)=X2+DX2
      P(2,8)=Y2-DY2
      P(3,8)=Z2
*
      GOTO 40
C
   30 CONTINUE
C
C             TRAP,PARA
C
C             Calculate the 8 vertex
C
      P(1,1)=-DZ*TX+TTH1*H1+TL1
      P(2,1)=+H1-DZ*TY
      P(3,1)=-DZ
      P(1,2)=-DZ*TX+TTH1*H1-TL1
      P(2,2)=+H1-DZ*TY
      P(3,2)=-DZ
      P(1,3)=-DZ*TX-TTH1*H1-BL1
      P(2,3)=-H1-DZ*TY
      P(3,3)=-DZ
      P(1,4)=-DZ*TX-TTH1*H1+BL1
      P(2,4)=-H1-DZ*TY
      P(3,4)=-DZ
      P(1,5)=+DZ*TX+TTH2*H2+TL2
      P(2,5)=+H2+DZ*TY
      P(3,5)=+DZ
      P(1,6)=+DZ*TX+TTH2*H2-TL2
      P(2,6)=+H2+DZ*TY
      P(3,6)=+DZ
      P(1,7)=+DZ*TX-TTH2*H2-BL2
      P(2,7)=-H2+DZ*TY
      P(3,7)=+DZ
      P(1,8)=+DZ*TX-TTH2*H2+BL2
      P(2,8)=-H2+DZ*TY
      P(3,8)=+DZ
C
   40 CONTINUE
C
C       BOX,TRD1,TRD2,TRAP,PARA --->>  call CGBOX
C
      IVCLOS=1
*SG
*  Size evaluation
*
      IF(IWORK.LE.2.AND.IVFUN.NE.0)THEN
*  NWB = n. words for each box
         NCGVOL=NCGVOL+NWB
         GOTO 999
      ENDIF
      ICPOIN=JCGOBJ+1
*   Creating object
*SG
      RMIN1=0
      RMIN2=0
      RMAX1=0
      RMAX2=0
      CALL CGBOX(P,4,4,300,Q(ICPOIN))
      DO 50  J=1,3
         T(J,1)=GRMAT(3*J-2,NLEVEL)
         T(J,2)=GRMAT(3*J-1,NLEVEL)
         T(J,3)=GRMAT(3*J,NLEVEL)
   50 CONTINUE
      CALL CGRIFL(T,Q(ICPOIN))
      CGERR=Q(ICPOIN)
      IF(CGERR.LE.0)THEN
         CALL GDCGER(CGERR)
         IF(KCGST.EQ.-2) GOTO 999
         IF(KCGST.EQ.-3) THEN
            KCGST=0
            WRITE(CHMAIL,10000)ISHAPE,NAMES(NLEVEL)
            CALL GMAIL(0,0)
            GOTO 999
         ENDIF
      ENDIF
      CALL CGCEV(1,Q(ICPOIN))
      DO 60  J=1,3
         T(J,1)=GRMAT(3*J-2,NLEVEL)
         T(J,2)=GRMAT(3*J-1,NLEVEL)
         T(J,3)=GRMAT(3*J,NLEVEL)
   60 CONTINUE
      CALL CGAFFI(T,Q(ICPOIN))
      XV=GTRAN(1,NLEVEL)
      YV=GTRAN(2,NLEVEL)
      ZV=GTRAN(3,NLEVEL)
      CALL CGSHIF(XV,YV,ZV,Q(ICPOIN))
***SG
*    Shifting object
      IF(KSHIFT.GT.0)THEN
         CALL GDSHIF(IVOLNA,ICPOIN)
      ENDIF
*
      IF(GBOOM.NE.0)THEN
         CALL GDBOMB(ICPOIN,ISHAPE)
         IF(ITSTCU.EQ.0)GOTO 999
      ENDIF
*
*
*
*   Hidden Volume Removal:
*   Computing volumes visibility and skipping
*   the unvisible ones; a great increase in speed
*   and a great reduction in n. of words used can be
*   obtained in this way.
*
      CALL CGMNMX(Q(ICPOIN),RRMIN,RRMAX)
      IF(ISUBLI.EQ.1)THEN
         AA1=RRMIN(1)-S1
         AA2=RRMIN(2)-S2
         AA3=RRMIN(3)-S3
         BB1=RRMAX(1)-SS1
         BB2=RRMAX(2)-SS2
         BB3=RRMAX(3)-SS3
         IF(AA1.GT.0.001.AND.AA2.GT.0.001.AND.AA3.GT.0.001.AND. BB1.LT.
     +   0.001.AND.BB2.LT.0.001.AND.BB3.LT.0.001)THEN
            IF(ISCOP.NE.1)THEN
               IF(IWORK.EQ.0.OR.IWORK.EQ.3)THEN
                  ITSTCU=0
                  NCGVOL=NCGVOL-NWB
                  GOTO 999
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(IPORLI.EQ.1)THEN
         S1=RRMIN(1)
         S2=RRMIN(2)
         S3=RRMIN(3)
         SS1=RRMAX(1)
         SS2=RRMAX(2)
         SS3=RRMAX(3)
         SRAGMX=0
         SRAGMN=0
         RAINT1=0
         RAINT2=0
      ENDIF
*   Create clipping objects
      IF(ICUT.NE.0.OR.IHOLE.NE.0)CALL GDCGSL(IVOLNA,ISHAPE)
*   Perspective view
      IF (IPRJ.EQ.IPERS) THEN
         CALL CGPERS(Q(ICPOIN))
      ENDIF
*   Inserting volumes in Hide + Wire Structures
      CALL GDCGHI(IWORK,NOBJ,NWWS,LSTEP,ISG,NTVOL,ISHAPE)
      GOTO 999
***SG
C
   70 CONTINUE
C
C     TUBE,CONE,TUBS,CONS  -----> call CGZREV
C
*
*     Checking Shape Parameters
*
      IF(RMIN1.GT.RMAX1) THEN
         WRITE(CHMAIL,10100)ISHAPE,NAMES(NLEVEL)
         CALL GMAIL(0,0)
      ENDIF
      IF(RMIN2.GT.RMAX2) THEN
         WRITE(CHMAIL,10200)ISHAPE,NAMES(NLEVEL)
         CALL GMAIL(0,0)
      ENDIF
      IF(PHIMIN.GT.PHIMAX)THEN
         WRITE(CHMAIL,10300)ISHAPE,NAMES(NLEVEL)
         CALL GMAIL(0,0)
      ENDIF
*
*
*   Checking if all Inner Radii are  0. ==> 'Closed' Volume
*
      IF(RMIN1.LE.0.00001.AND.RMIN2.LE.0.00001)IVCLOS=1
*SG
*   Size evaluation
      IF(IWORK.LE.2.AND.IVFUN.NE.0)THEN
*   NWPROD = n. words for each body of revolution
         NWPROD=NWREV*(NANG+1)
         NCGVOL=NCGVOL+NWPROD
         GOTO 999
      ENDIF
*  Creating object
      ICPOIN=JCGOBJ+1
*SG
      IF(ISHAPE.EQ.29)THEN
         SAL=PAR(8)
         IF(PAR(11).GT.SAL)SAL=PAR(11)
         PAR3=MAX(PAR(3),0.)
         Z2=PAR3+1.001*RMAX1*SQRT((1-SAL*SAL)/(SAL*SAL))
      ENDIF
      XZ(1,1)=RMIN1
      XZ(2,1)=-Z2
      XZ(1,2)=RMAX1
      XZ(2,2)=-Z2
      XZ(1,3)=RMAX2
      XZ(2,3)=Z2
      XZ(1,4)=RMIN2
      XZ(2,4)=Z2
      CALL CGZREV(XZ,ANG1,ANG2,NANG,4000,Q(ICPOIN))
      DO 80  J=1,3
         T(J,1)=GRMAT(3*J-2,NLEVEL)
         T(J,2)=GRMAT(3*J-1,NLEVEL)
         T(J,3)=GRMAT(3*J,NLEVEL)
   80 CONTINUE
      CALL CGRIFL(T,Q(ICPOIN))
      CGERR=Q(ICPOIN)
      IF(CGERR.LE.0)THEN
         CALL GDCGER(CGERR)
         IF(KCGST.EQ.-2) GOTO 999
         IF(KCGST.EQ.-3) THEN
            KCGST=0
            WRITE(CHMAIL,10000)ISHAPE,NAMES(NLEVEL)
            CALL GMAIL(0,0)
            GOTO 999
         ENDIF
      ENDIF
      CALL CGCEV(1,Q(ICPOIN))
      IF(ISHAPE.EQ.29)THEN
         SLI1(1)=-PAR(6)
         SLI1(2)=-PAR(7)
         SLI1(3)=-PAR(8)
         SLI1(4)=-PAR(3)*PAR(8)
         SLI2(1)=-PAR(9)
         SLI2(2)=-PAR(10)
         SLI2(3)=-PAR(11)
         SLI2(4)=+PAR(3)*PAR(11)
         ISL1=JCGOBJ+4000
         CALL CGSLIC(Q(ICPOIN),SLI1,4000,Q(ISL1))
         ISL2=JCGOBJ+8000
         CALL CGSLIC(Q(ISL1),SLI2,4000,Q(ISL2))
         ICPOIN=ISL2
         CALL CGCEV(1,Q(ICPOIN))
      ENDIF
      DO 90  J=1,3
         T(J,1)=GRMAT(3*J-2,NLEVEL)
         T(J,2)=GRMAT(3*J-1,NLEVEL)
         T(J,3)=GRMAT(3*J,NLEVEL)
   90 CONTINUE
      CALL CGAFFI(T,Q(ICPOIN))
      XV=GTRAN(1,NLEVEL)
      YV=GTRAN(2,NLEVEL)
      ZV=GTRAN(3,NLEVEL)
      CALL CGSHIF(XV,YV,ZV,Q(ICPOIN))
***SG
*    Shifting object
      IF(KSHIFT.GT.0)THEN
         CALL GDSHIF(IVOLNA,ICPOIN)
      ENDIF
*
      IF(GBOOM.NE.0)THEN
         CALL GDBOMB(ICPOIN,ISHAPE)
         IF(ITSTCU.EQ.0)GOTO 999
      ENDIF
*
*
*   Hidden Volume Removal:
*   Computing closed volumes visibility and skipping
*   the unvisible ones; a great increase in speed
*   and a great reduction in n. of words used are obtained
*   in this way.
*
      CALL CGMNMX(Q(ICPOIN),RRMIN,RRMAX)
      IF(ISUBLI.EQ.1)THEN
         AA1=RRMIN(1)-S1
         AA2=RRMIN(2)-S2
         AA3=RRMIN(3)-S3
         BB1=RRMAX(1)-SS1
         BB2=RRMAX(2)-SS2
         BB3=RRMAX(3)-SS3
         IF(AA1.GT.0.001.AND.AA2.GT.0.001.AND.AA3.GT.0.001.AND. BB1.LT.
     +   -0.001.AND.BB2.LT.-0.001.AND.BB3.LT.-0.001)THEN
            IF(ISHAPE.EQ.7.OR.ISHAPE.EQ.8)THEN
               IF((RMAX2.LT.SRAGMX.AND.RMAX1.LT.SRAGMN).OR. (SRAGMX.EQ.
     +         0))THEN
                  IF((RMIN2.GT.RAINT2.AND.RMIN1.GT.RAINT1).OR. (RAINT2.
     +            EQ.0))THEN
                     IF(IWORK.EQ.0.OR.IWORK.EQ.3)THEN
                        ITSTCU=0
                        NCGVOL=NCGVOL-NWPROD
                        GOTO 999
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF(SRAGMX.NE.0.)THEN
               DO 100 ITER=1,IPORNT
                  IF(RMAX1.EQ.PORMAR(ITER))GOTO 110
                  IF(RMIN1.EQ.PORMIR(ITER))THEN
                     IF(PORMIR(ITER).NE.0.)GOTO 110
                  ENDIF
  100          CONTINUE
            ENDIF
            IF(IWORK.EQ.0.OR.IWORK.EQ.3)THEN
               ITSTCU=0
               NCGVOL=NCGVOL-NWPROD
               GOTO 999
            ENDIF
         ENDIF
      ENDIF
      IF(IPORLI.EQ.1)THEN
         S1=RRMIN(1)
         S2=RRMIN(2)
         S3=RRMIN(3)
         SS1=RRMAX(1)
         SS2=RRMAX(2)
         SS3=RRMAX(3)
         SRAGMX=RMAX2
         SRAGMN=RMAX1
         RAINT1=RMIN1
         RAINT2=RMIN2
         IPORNT=1
         PORMAR(1)=RMAX2
         PORMIR(1)=RMIN1
      ENDIF
  110 CONTINUE
*   Create clipping objects
      IF(ICUT.NE.0.OR.IHOLE.NE.0)CALL GDCGSL(IVOLNA,ISHAPE)
*   Perspective view
      IF (IPRJ.EQ.IPERS) THEN
         CALL CGPERS(Q(ICPOIN))
      ENDIF
*   Inserting objects in Hide + Wire structures
      CALL GDCGHI(IWORK,NOBJ,NWWS,LSTEP,ISG,NTVOL,ISHAPE)
      GOTO 999
***SG
C
  120 CONTINUE
C
C     SPHE  -----> call CGSPHE
C
      IVCLOS=1
*SG
*   Size evaluation
      IF(IWORK.LE.2.AND.IVFUN.NE.0)THEN
*  NWS = n. words for each sphere
         NCGVOL=NCGVOL+NWS
         GOTO 999
      ENDIF
*
      R=RMAX
      RMAX2=R
      RMAX1=0
      RMIN1=0
      RMIN2=0
      NLAT=11
      NLON=11
      NWOR=4000
      IF(IWORK.EQ.3.AND.(PAR(3).EQ.0.AND.(PAR(4).EQ.0.OR.
     +PAR(4).EQ.180)))THEN
         NLAT=29
         NLON=29
         NWOR=30000
      ENDIF
      ICPOIN=JCGOBJ+1
*  Creating object
      CALL CGSPHE(R,NLAT,NLON,NWOR,Q(ICPOIN))
      DO 130 J=1,3
         T(J,1)=GRMAT(3*J-2,NLEVEL)
         T(J,2)=GRMAT(3*J-1,NLEVEL)
         T(J,3)=GRMAT(3*J,NLEVEL)
  130 CONTINUE
      CALL CGRIFL(T,Q(ICPOIN))
*SG
      CGERR=Q(ICPOIN)
      IF(CGERR.LE.0)THEN
         CALL GDCGER(CGERR)
         IF(KCGST.EQ.-2) GOTO 999
         IF(KCGST.EQ.-3) THEN
            KCGST=0
            WRITE(CHMAIL,10000)ISHAPE,NAMES(NLEVEL)
            CALL GMAIL(0,0)
            GOTO 999
         ENDIF
      ENDIF
      CALL CGCEV(1,Q(ICPOIN))
      IF(PAR(3).NE.0.OR.(PAR(4).NE.0.AND.PAR(4).NE.180))THEN
         ISHAPE=99
         SPI1(1)=-COS((90-PAR(3))*DEGRAD)
         SPI1(2)=0
         SPI1(3)=-COS(PAR(3)*DEGRAD)
         SPI1(4)=0
         SPI2(1)=-COS((90-PAR(4))*DEGRAD)
         SPI2(2)=0
         SPI2(3)=-COS(PAR(4)*DEGRAD)
         SPI2(4)=0
         ISP1=JCGOBJ+4000
         CALL CGSLIC(Q(ICPOIN),SPI1,4000,Q(ISP1))
         ISP2=JCGOBJ+8000
         CALL CGSLIC(Q(ISP1),SPI2,4000,Q(ISP2))
         ICPOIN=ISP2
         CALL CGCEV(1,Q(ICPOIN))
      ENDIF
      DO 140 J=1,3
         T(J,1)=GRMAT(3*J-2,NLEVEL)
         T(J,2)=GRMAT(3*J-1,NLEVEL)
         T(J,3)=GRMAT(3*J,NLEVEL)
  140 CONTINUE
      CALL CGAFFI(T,Q(ICPOIN))
      XV=GTRAN(1,NLEVEL)
      YV=GTRAN(2,NLEVEL)
      ZV=GTRAN(3,NLEVEL)
      CALL CGSHIF(XV,YV,ZV,Q(ICPOIN))
***SG
*    Shifting object
      IF(KSHIFT.GT.0)THEN
         CALL GDSHIF(IVOLNA,ICPOIN)
      ENDIF
*
      IF(GBOOM.NE.0)THEN
         CALL GDBOMB(ICPOIN,ISHAPE)
         IF(ITSTCU.EQ.0)GOTO 999
      ENDIF
*
*
*   Hidden Volume Removal:
*   Computing closed volumes visibility and skipping
*   the unvisible ones; a great increase in speed
*   and a great reduction in n. of words used are obtained
*   in this way.
*
      CALL CGMNMX(Q(ICPOIN),RRMIN,RRMAX)
      IF(ISUBLI.EQ.1)THEN
         AA1=RRMIN(1)-S1
         AA2=RRMIN(2)-S2
         AA3=RRMIN(3)-S3
         BB1=RRMAX(1)-SS1
         BB2=RRMAX(2)-SS2
         BB3=RRMAX(3)-SS3
         IF(AA1.GT.0.001.AND.AA2.GT.0.001.AND.AA3.GT.0.001.AND. BB1.LT.
     +   -0.001.AND.BB2.LT.-0.001.AND.BB3.LT.-0.001)THEN
            IF(ISHAPE.NE.99)THEN
               IF(IWORK.EQ.0.OR.IWORK.EQ.3)THEN
                  ITSTCU=0
                  NCGVOL=NCGVOL-NWS
                  GOTO 999
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(IPORLI.EQ.1)THEN
         S1=RRMIN(1)
         S2=RRMIN(2)
         S3=RRMIN(3)
         SS1=RRMAX(1)
         SS2=RRMAX(2)
         SS3=RRMAX(3)
         SRAGMX=R
         SRAGMN=0.
         RAINT1=0.
         RAINT2=0.
         IPORNT=1
         PORMAR(1)=R
         PORMIR(1)=0.
      ENDIF
*   Create clipping objects
      IF(ICUT.NE.0.OR.IHOLE.NE.0)CALL GDCGSL(IVOLNA,ISHAPE)
*   Perspective view
      IF (IPRJ.EQ.IPERS) THEN
         CALL CGPERS(Q(ICPOIN))
      ENDIF
*   Inserting objects in Hide + Wire structures
      CALL GDCGHI(IWORK,NOBJ,NWWS,LSTEP,ISG,NTVOL,ISHAPE)
      GOTO 999
***SG
*
  150 CONTINUE
C
C         PGON   ---->  call CGZREV
C
      NTVOL=NZ-1
      ANG1=PHIMIN
      ANG2=PHIMAX
**SG
      AANG=ABS(ANG2-ANG1)
      AZLAT=AANG*APPROS
      ZLAT=AZLAT/360
      NANG=ZLAT
      IF(NANG.EQ.0)NANG=1
      AZ=ZLAT-NANG
      IF(AZ.GT..5)NANG=NANG+1
      IF(NDIVAN.LT.NANG)THEN
         NANG=NDIVAN
*         WRITE(CHMAIL,10400)NANG,NAMES(NLEVEL)
*         CALL GMAIL(0,0)
      ENDIF
      AATMAX=NANG*360./AANG
      LATMAX=AATMAX
      ALA=AATMAX-LATMAX
      IF(ALA.GT..5)LATMAX=LATMAX+1
**SG
      AFINV=1./COS(PI/LATMAX)
      FINV=ABS(AFINV)
      JSURZ=1
      ZR(1)=PAR(5)
      RMIR(1)=PAR(6)*FINV
      RMAR(1)=PAR(7)*FINV
*SG
      RMAR(1)=RMAR(1)+.001
*SG
      DO 160 I=1,NTVOL
*         ZA=PAR(5+3*(I-1))
         ZB=PAR(5+3*I)
**SG
         ZB=ZB+.001
*********         DIFZ=ABS(ZB-ZA)
*********         IF(DIFZ.LT.0.001)GOTO 220
**SG
         JSURZ=JSURZ+1
         ZR(JSURZ)=ZB
         RMIR(JSURZ)=PAR(6+3*I)*FINV
         RMAR(JSURZ)=PAR(7+3*I)*FINV
**SG
         RMAR(JSURZ)=RMAR(JSURZ)+.001
*
  160 CONTINUE
*
*   Checking if all Inner Radii are  0. ==> 'Closed' Volume
*
*      NRAD=NTVOL+1
*      DO 230 I=1,NRAD
*         IF(RMIR(I).GT.0.00001)GOTO 240
*  230 CONTINUE
*      IVCLOS=1
*  240 CONTINUE
*
*   Size evaluation
      IF(IWORK.LE.2.AND.IVFUN.NE.0)THEN
         NCGVOL=NCGVOL+NWREV*(NANG+1)*NTVOL
         GOTO 999
      ENDIF
      IF(IPORLI.EQ.1)THEN
         SRAGMN=10000.
         RAINT1=10000.
      ENDIF
*
      DO 220  IVOL=1,NTVOL
         ISG=ISG+1
         IVCLOS=1
         IF((RMIR(IVOL).GT.0.00001).OR.(RMIR(IVOL+1).GT.
     +   0.00001))IVCLOS=0
         XZ(1,1)=RMIR(IVOL)
         XZ(2,1)=ZR(IVOL)
         XZ(1,2)=RMAR(IVOL)
         XZ(2,2)=ZR(IVOL)
         XZ(1,3)=RMAR(IVOL+1)
         XZ(2,3)=ZR(IVOL+1)
         XZ(1,4)=RMIR(IVOL+1)
         XZ(2,4)=ZR(IVOL+1)
         ZR(IVOL+1)=ZR(IVOL+1)+.001
         ICPOIN=JCGOBJ+1
*   Creating object
**SG
         CALL CGZREV(XZ,ANG1,ANG2,NANG,4000,Q(ICPOIN))
         DO 170 J=1,3
            T(J,1)=GRMAT(3*J-2,NLEVEL)
            T(J,2)=GRMAT(3*J-1,NLEVEL)
            T(J,3)=GRMAT(3*J,NLEVEL)
  170    CONTINUE
         CALL CGRIFL(T,Q(ICPOIN))
         CGERR=Q(ICPOIN)
         IF(CGERR.LE.0)THEN
            CALL GDCGER(CGERR)
            IF(KCGST.EQ.-2) GOTO 999
            IF(KCGST.EQ.-3) THEN
               KCGST=0
               WRITE(CHMAIL,10000)ISHAPE,NAMES(NLEVEL)
               CALL GMAIL(0,0)
               WRITE(CHMAIL,10400)(PAR(I),I=1,4)
               CALL GMAIL(0,0)
               DO 180 J=1,NZ
                  ZPR=PAR(5+(J-1)*3)
                  RMIPR=PAR(6+(J-1)*3)
                  RMAPR=PAR(7+(J-1)*3)
                  WRITE(CHMAIL,10600)J,ZPR,RMIPR,RMAPR
                  CALL GMAIL(0,0)
  180          CONTINUE
               GOTO 999
            ENDIF
         ENDIF
         CALL CGCEV(1,Q(ICPOIN))
         DO 190  J=1,3
            T(J,1)=GRMAT(3*J-2,NLEVEL)
            T(J,2)=GRMAT(3*J-1,NLEVEL)
            T(J,3)=GRMAT(3*J,NLEVEL)
  190    CONTINUE
         CALL CGAFFI(T,Q(ICPOIN))
         XV=GTRAN(1,NLEVEL)
         YV=GTRAN(2,NLEVEL)
         ZV=GTRAN(3,NLEVEL)
         CALL CGSHIF(XV,YV,ZV,Q(ICPOIN))
***SG
*    Shifting object
         IF(KSHIFT.GT.0)THEN
            CALL GDSHIF(IVOLNA,ICPOIN)
         ENDIF
*
         IF(GBOOM.NE.0)THEN
            CALL GDBOMB(ICPOIN,ISHAPE)
            IF(ITSTCU.EQ.0)GOTO 220
         ENDIF
*
*
*   Hidden Volume Removal:
*   Computing closed volumes visibility and skipping
*   the unvisible ones; a great increase in speed
*   and a great reduction in n. of words used are obtained
*   in this way.
*
         CALL CGMNMX(Q(ICPOIN),RRMIN,RRMAX)
         IF(ISUBLI.EQ.1)THEN
            AA1=RRMIN(1)-S1
            AA2=RRMIN(2)-S2
            AA3=RRMIN(3)-S3
            BB1=RRMAX(1)-SS1
            BB2=RRMAX(2)-SS2
            BB3=RRMAX(3)-SS3
            IF(AA1.GT.0.001.AND.AA2.GT.0.001.AND.AA3.GT.0.001.AND.
     +      BB1.LT.-0.001.AND.BB2.LT.-0.001.AND.BB3.LT.-0.001)THEN
               AMARMA(IVOL) =MIN(RMAR(IVOL),RMAR(IVOL+1))
               AMARMA(IVOL+1)=MAX(RMAR(IVOL),RMAR(IVOL+1))
               AMIRMA(IVOL) =MIN(RMIR(IVOL),RMIR(IVOL+1))
               AMIRMA(IVOL+1)=MAX(RMIR(IVOL),RMIR(IVOL+1))
               RMAX1=AMARMA(IVOL)
               RMAX2=AMARMA(IVOL+1)
               RMIN1=AMIRMA(IVOL)
               RMIN2=AMIRMA(IVOL+1)
               IF(SRAGMX.NE.0.)THEN
                  DO 200 ITER=1,IPORNT
                     IF(RMAX1.EQ.PORMAR(ITER).OR.RMAX2.EQ.PORMAR(ITER))
     +               GOTO 210
                     IF(RMIN1.EQ.PORMIR(ITER).OR.RMIN2.EQ.PORMIR(ITER))
     +               THEN
                        IF(PORMIR(ITER).NE.0.)GOTO 210
                     ENDIF
  200             CONTINUE
               ENDIF
               IF(ISCOP.EQ.1)THEN
                  IF((AMARMA(IVOL+1).GE.SRAGMX.OR.AMARMA(IVOL) .GT.SRAG
     +            MN))GOTO 210
                  IF((AMIRMA(IVOL+1).LE.RAINT2.OR.AMIRMA(IVOL)
     +            .LE.RAINT1).AND.(RAINT2.NE.0))GOTO 210
               ENDIF
               IF(IWORK.EQ.0.OR.IWORK.EQ.3)THEN
                  ITSTCU=0
                  NCGVOL=NCGVOL-NWPROD
                  GOTO 220
               ENDIF
            ENDIF
         ENDIF
         IF(IPORLI.EQ.1)THEN
            IF(RRMIN(1).LT.S1)S1=RRMIN(1)
            IF(RRMIN(2).LT.S2)S2=RRMIN(2)
            IF(RRMIN(3).LT.S3)S3=RRMIN(3)
            IF(RRMAX(1).GT.SS1)SS1=RRMAX(1)
            IF(RRMAX(2).GT.SS2)SS2=RRMAX(2)
            IF(RRMAX(3).GT.SS3)SS3=RRMAX(3)
            IF(RMAR(IVOL).GT.SRAGMX)SRAGMX=RMAR(IVOL)
            IF(RMAR(IVOL).LT.SRAGMN)SRAGMN=RMAR(IVOL)
            IF(RMAR(IVOL+1).GT.SRAGMX)SRAGMX=RMAR(IVOL+1)
            IF(RMAR(IVOL+1).LT.SRAGMN)SRAGMN=RMAR(IVOL+1)
            IF(RMIR(IVOL).GT.RAINT2)RAINT2=RMIR(IVOL)
            IF(RMIR(IVOL).LT.RAINT1)RAINT1=RMIR(IVOL)
            IF(RMIR(IVOL+1).GT.RAINT2)RAINT2=RMIR(IVOL+1)
            IF(RMIR(IVOL+1).LT.RAINT1)RAINT1=RMIR(IVOL+1)
            PORMAR(IVOL)=RMAR(IVOL)
            PORMIR(IVOL)=RMIR(IVOL)
            IPORNT =NTVOL
         ENDIF
  210    CONTINUE
*   Create clipping objects
         IF(ICUT.NE.0.OR.IHOLE.NE.0)CALL GDCGSL(IVOLNA,ISHAPE)
*   Perspective view
         IF (IPRJ.EQ.IPERS) THEN
            CALL CGPERS(Q(ICPOIN))
         ENDIF
*   Inserting objects in Hide + Wire structures
         CALL GDCGHI(IWORK,NOBJ,NWWS,LSTEP,ISG,NTVOL,ISHAPE)
***SG
  220 CONTINUE
      GOTO 999
*
  230 CONTINUE
C
C         PCON    ---->  call CGZREV
C
      NTVOL=NZ-1
      ANG1=PHIMIN
      ANG2=PHIMAX
**SG
      AANG=ABS(ANG2-ANG1)
      AZLAT=AANG*APPROS
      ZLAT=AZLAT/360
      NANG=ZLAT
      IF(NANG.EQ.0)NANG=1
      AZ=ZLAT-NANG
      IF(AZ.GT..5)NANG=NANG+1
**SG
      AFINV=1./COS(PI/APPROS)
      FINV=ABS(AFINV)
      JSURZ=1
      ZR(1)=PAR(4)
      RMIR(1)=PAR(5)*FINV
      RMAR(1)=PAR(6)*FINV
*SG
      RMAR(1)=RMAR(1)+.1
*SG
      DO 240 I=1,NTVOL
*         ZA=PAR(4+3*(I-1))
         ZB=PAR(4+3*I)
**SG
         ZB=ZB+.001
********         DIFZ=ABS(ZB-ZA)
********         IF(DIFZ.LT.0.001)GOTO 290
**SG
         JSURZ=JSURZ+1
         ZR(JSURZ)=ZB
         RMIR(JSURZ)=PAR(5+3*I)*FINV
         RMAR(JSURZ)=PAR(6+3*I)*FINV
**SG
         RMAR(JSURZ)=RMAR(JSURZ)+.1
*
  240 CONTINUE
*
*   Checking if all Inner Radii are  0. ==> 'Closed' Volume
*
*      NRAD=NTVOL+1
*      DO 300 I=1,NRAD
*         IF(RMIR(I).GT.0.00001)GOTO 310
*  300 CONTINUE
*      IVCLOS=1
*  310 CONTINUE
*
*
*   Size evaluation
      IF(IWORK.LE.2.AND.IVFUN.NE.0)THEN
         NCGVOL=NCGVOL+NWREV*(NANG+1)*NTVOL
         GOTO 999
      ENDIF
      IF(IPORLI.EQ.1)THEN
         SRAGMN=10000.
         RAINT1=10000.
      ENDIF
*
      DO 300 IVOL=1,NTVOL
         ISG=ISG+1
         IVCLOS=1
         IF((RMIR(IVOL).GT.0.00001).OR.(RMIR(IVOL+1).GT.
     +   0.00001))IVCLOS=0
         XZ(1,1)=RMIR(IVOL)
         XZ(2,1)=ZR(IVOL)
         XZ(1,2)=RMAR(IVOL)
         XZ(2,2)=ZR(IVOL)
         XZ(1,3)=RMAR(IVOL+1)
         XZ(2,3)=ZR(IVOL+1)
         XZ(1,4)=RMIR(IVOL+1)
         XZ(2,4)=ZR(IVOL+1)
         ZR(IVOL+1)=ZR(IVOL+1)+.1
         ICPOIN=JCGOBJ+1
*  Creating object
**SG
         CALL CGZREV(XZ,ANG1,ANG2,NANG,4000,Q(ICPOIN))
         DO 250 J=1,3
            T(J,1)=GRMAT(3*J-2,NLEVEL)
            T(J,2)=GRMAT(3*J-1,NLEVEL)
            T(J,3)=GRMAT(3*J,NLEVEL)
  250    CONTINUE
         CALL CGRIFL(T,Q(ICPOIN))
         CALL CGCEV(-1,Q(ICPOIN))
         CGERR=Q(ICPOIN)
         IF(CGERR.LE.0)THEN
            CALL GDCGER(CGERR)
            IF(KCGST.EQ.-2) GOTO 999
            IF(KCGST.EQ.-3) THEN
               KCGST=0
               WRITE(CHMAIL,10000)ISHAPE,NAMES(NLEVEL)
               CALL GMAIL(0,0)
               WRITE(CHMAIL,10500)(PAR(I),I=1,3)
               CALL GMAIL(0,0)
               DO 260 J=1,NZ
                  ZPR=PAR(4+(J-1)*3)
                  RMIPR=PAR(5+(J-1)*3)
                  RMAPR=PAR(6+(J-1)*3)
                  WRITE(CHMAIL,10600)J,ZPR,RMIPR,RMAPR
                  CALL GMAIL(0,0)
  260          CONTINUE
               GOTO 999
            ENDIF
         ENDIF
         CALL CGCEV(1,Q(ICPOIN))
         DO 270 J=1,3
            T(1,J)=GRMAT(3*J-2,NLEVEL)
            T(2,J)=GRMAT(3*J-1,NLEVEL)
            T(3,J)=GRMAT(3*J,NLEVEL)
            T(4,J)=0.
  270    CONTINUE
         CALL CGAFFI(T,Q(ICPOIN))
         XV=GTRAN(1,NLEVEL)
         YV=GTRAN(2,NLEVEL)
         ZV=GTRAN(3,NLEVEL)
         CALL CGSHIF(XV,YV,ZV,Q(ICPOIN))
*
***SG
*    Shifting object
         IF(KSHIFT.GT.0)THEN
            CALL GDSHIF(IVOLNA,ICPOIN)
         ENDIF
*
         IF(GBOOM.NE.0)THEN
            CALL GDBOMB(ICPOIN,ISHAPE)
            IF(ITSTCU.EQ.0)GOTO 300
         ENDIF
*
*
*   Hidden Volume Removal:
*   Computing closed volumes visibility and skipping
*   the unvisible ones; a great increase in speed
*   and a great reduction in n. of words used are obtained
*   in this way.
*
         CALL CGMNMX(Q(ICPOIN),RRMIN,RRMAX)
         IF(ISUBLI.EQ.1)THEN
            AA1=RRMIN(1)-S1
            AA2=RRMIN(2)-S2
            AA3=RRMIN(3)-S3
            BB1=RRMAX(1)-SS1
            BB2=RRMAX(2)-SS2
            BB3=RRMAX(3)-SS3
            IF(AA1.GT.0.001.AND.AA2.GT.0.001.AND.AA3.GT.0.001.AND.
     +      BB1.LT.-0.001.AND.BB2.LT.-0.001.AND.BB3.LT.-0.001)THEN
               AMARMA(IVOL) =MIN(RMAR(IVOL),RMAR(IVOL+1))
               AMARMA(IVOL+1)=MAX(RMAR(IVOL),RMAR(IVOL+1))
               AMIRMA(IVOL) =MIN(RMIR(IVOL),RMIR(IVOL+1))
               AMIRMA(IVOL+1)=MAX(RMIR(IVOL),RMIR(IVOL+1))
               RMAX1=AMARMA(IVOL)
               RMAX2=AMARMA(IVOL+1)
               RMIN1=AMIRMA(IVOL)
               RMIN2=AMIRMA(IVOL+1)
               IF(SRAGMX.NE.0.)THEN
                  DO 280 ITER=1,IPORNT
                     IF(RMAX1.EQ.PORMAR(ITER).OR.RMAX2.EQ.PORMAR(ITER))
     +               GOTO 290
                     IF(RMIN1.EQ.PORMIR(ITER).OR.RMIN2.EQ.PORMIR(ITER))
     +               THEN
                        IF(PORMIR(ITER).NE.0)GOTO 290
                     ENDIF
  280             CONTINUE
               ENDIF
               IF(ISCOP.EQ.1)THEN
                  IF((AMARMA(IVOL+1).GE.SRAGMX.OR.AMARMA(IVOL) .GE.SRAG
     +            MN))GOTO 290
                  IF((AMIRMA(IVOL+1).LE.RAINT2.AND.AMIRMA(IVOL)
     +            .LE.RAINT1).AND.(RAINT2.NE.0))GOTO 290
               ENDIF
               IF(IWORK.EQ.0.OR.IWORK.EQ.3)THEN
                  NCGVOL=NCGVOL-NWPROD
                  ITSTCU=0
                  GOTO 300
               ENDIF
            ENDIF
         ENDIF
         IF(IPORLI.EQ.1)THEN
            IF(RRMIN(1).LT.S1)S1=RRMIN(1)
            IF(RRMIN(2).LT.S2)S2=RRMIN(2)
            IF(RRMIN(3).LT.S3)S3=RRMIN(3)
            IF(RRMAX(1).GT.SS1)SS1=RRMAX(1)
            IF(RRMAX(2).GT.SS2)SS2=RRMAX(2)
            IF(RRMAX(3).GT.SS3)SS3=RRMAX(3)
            IF(RMAR(IVOL).GT.SRAGMX)SRAGMX=RMAR(IVOL)
            IF(RMAR(IVOL).LT.SRAGMN)SRAGMN=RMAR(IVOL)
            IF(RMAR(IVOL+1).GT.SRAGMX)SRAGMX=RMAR(IVOL+1)
            IF(RMAR(IVOL+1).LT.SRAGMN)SRAGMN=RMAR(IVOL+1)
            IF(RMIR(IVOL).GT.RAINT2)RAINT2=RMIR(IVOL)
            IF(RMIR(IVOL).LT.RAINT1)RAINT1=RMIR(IVOL)
            IF(RMIR(IVOL+1).GT.RAINT2)RAINT2=RMIR(IVOL+1)
            IF(RMIR(IVOL+1).LT.RAINT1)RAINT1=RMIR(IVOL+1)
            PORMAR(IVOL)=RMAR(IVOL)
            PORMIR(IVOL)=RMIR(IVOL)
            IPORNT =NTVOL
         ENDIF
  290    CONTINUE
*   Create clipping objects
         IF(ICUT.NE.0.OR.IHOLE.NE.0)CALL GDCGSL(IVOLNA,ISHAPE)
*   Perspective view
         IF (IPRJ.EQ.IPERS) THEN
            CALL CGPERS(Q(ICPOIN))
         ENDIF
*   Inserting object in Hide + Wire structures
         CALL GDCGHI(IWORK,NOBJ,NWWS,LSTEP,ISG,NTVOL,ISHAPE)
  300 CONTINUE
      GOTO 999
*
***SG
*
10000   FORMAT(' Check Parameters of Shape ',I3,' in volume ',A4)
10100   FORMAT(' Warning >>> RMIN1 greater than RMAX1 for shape '
     +  ,I3,' in volume ',A4)
10200   FORMAT(' Warning >>> RMIN2 greater than RMAX2 for shape '
     +  ,I3,' in volume ',A4)
10300   FORMAT(' Warning >>> PHIMIN greater than PHIMAX for shape'
     +  ,I3,' in volume ',A4)
*10400   FORMAT(' PGON with NPDV = ',I5,' in volume ',A4,' NPDV very
*     +  large . It must be < 30 . Volume will not be drawn. ')
10400   FORMAT(' PHI1 = ',F8.3,' DPHI = ',F8.3,' NPDV = ',F8.1
     +  ,' NZ   = ',F8.1)
10500   FORMAT(' PHI1 = ',F8.3,' DPHI = ',F8.3,' NZ  = ',F8.1)
10600   FORMAT(' J = ',I5,' Z = ',F8.3,' RMIN = ',F8.3
     +  ,' RMAX = ',F8.3)
*10800   FORMAT(' Please, increase size of Zebra store by ',I10,
*     +         ' words')
*
***SG
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:03  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 10/01/95  16.43.07  by  S.Giani
*-- Author :
      SUBROUTINE GDPRTR(NAME,YLEV,XLEV,SIZE,FACHV,
     +IDUP,IDRAW,IONL,IDET,ISEEN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draws a box centered in (XLEV,YLEV) of size SIZE.        *
C.    *       On the top half box draws its name NAME and, if IDUP<>1, *
C.    *       on the bottom half draw IDUP preceded by a symbol        *
C.    *       (IDUP<0 means division, IDUP>1 means multiple position,  *
C.    *       IDUP=1 means single position).                           *
C.    *       If IDRAW=1 then over the box draws a normal-view of      *
C.    *       volume NAME (with same seen attributes set by GDSPEC).   *
C.    *       If IONL<>1 then draw two other boxes, slightly shifted,  *
C.    *       behind the original one.                                 *
C.    *       If IDET=1 then draw another box, slightly larger,        *
C.    *       around the original one.                                 *
C.    *       If ISEEN=0 then draw the box(es) with dashed lines.      *
C.    *                                                                *
C.    *    ==>Called by : GDTREE                                       *
C.    *       Author : P.Zanarini ; S.Giani  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcdraw.inc"
      DIMENSION X(5),Y(5)
      DIMENSION NAMT0(16),NAM0(4)
      DIMENSION PAR(100)
      DIMENSION NNDM(14),INDM(4,14)
      DIMENSION ATT(10)
      CHARACTER*4 NAME,NAMSEE,ISON
      CHARACTER*14 CHNAM0
      SAVE NNDM,INDM
      COMMON/NPILEV/NPILE
C
      DATA NNDM/0,0,0,4,0,2,0,2,4,3,4,3,0,0/
      DATA INDM/0,0,0,0, 0,0,0,0, 0,0,0,0, 2,3,7,11, 0,0,0,0,
     +          4,5,0,0, 0,0,0,0, 6,7,0,0, 3,4,5,6, 4,5,6,0,
     +          1,2,3,4, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C.
C.    ------------------------------------------------------------------
C.
      IF (ISEEN.EQ.0) THEN
         LINSTY=4
      ELSE
         LINSTY=1
      ENDIF
      CALL MVBITS(LINSTY,0,3,LINATT,10)
C
C            Compute box vertex
C
      SIZ2=SIZE/2.
      SIZ4=SIZE/4.
      SIZ8=SIZE/8.
***      X(1)=XLEV-SIZ2
      X(1)=XLEV-(SIZ2*2.*FACHV)
      X(2)=X(1)
***      X(3)=XLEV+SIZ2
      X(3)=XLEV+(SIZ2*2.*FACHV)
      X(4)=X(3)
      X(5)=X(1)
      Y(1)=YLEV+SIZ2
      Y(2)=YLEV-SIZ2
      Y(3)=Y(2)
      Y(4)=Y(1)
      Y(5)=Y(1)
C
C            Draw box
C
      IWTY=IGIWTY(1)
      CALL IGQWK(0,'HWCO',R)
      IF(IWTY.GT.10.OR.IWTY.LT.1.OR.R.EQ.2.)THEN
         LINCOL=1
         CALL MVBITS(LINCOL,0,8,LINATT,16)
         CALL GDRAWV(X,Y,5)
      ELSE
         CALL IGSET('LWID',1.)
         CALL UCTOH(NAME,IIIN,4,4)
         CALL IGPID(1,'Tree',IIIN,' ')
         CALL IGPAVE(X(1),X(3),Y(1),Y(3),SIZ8,1004,1001,'TRS')
      ENDIF
C
C            If volume is also detector draw another box outside
C
      IF (IDET.EQ.1) THEN
         X(1)=X(1)-SIZ8
         X(2)=X(2)-SIZ8
         X(3)=X(3)+SIZ8
         X(4)=X(4)+SIZ8
         X(5)=X(5)-SIZ8
         Y(1)=Y(1)+SIZ8
         Y(2)=Y(2)-SIZ8
         Y(3)=Y(3)-SIZ8
         Y(4)=Y(4)+SIZ8
         Y(5)=Y(5)+SIZ8
         LINCOL=1
         CALL MVBITS(LINCOL,0,8,LINATT,16)
         CALL GDRAWV(X,Y,5)
      ENDIF
C
C            If volume is not 'ONLY' draw two other boxes behind
C
      IF (IONL.NE.1) THEN
         X(1)=XLEV+(SIZ2*2.*FACHV)-SIZ8
         Y(1)=YLEV+SIZ2
         X(2)=X(1)
         Y(2)=Y(1)+SIZ8
         X(3)=X(2)-(SIZE*2.*FACHV)
         Y(3)=Y(2)
         X(4)=X(3)
         Y(4)=Y(3)-SIZE
         X(5)=X(4)+SIZ8
         Y(5)=Y(4)
         LINCOL=1
         CALL MVBITS(LINCOL,0,8,LINATT,16)
         CALL GDRAWV(X,Y,5)
         DO 10 I=1,5
            X(I)=X(I)-SIZ8
            Y(I)=Y(I)+SIZ8
   10    CONTINUE
         CALL GDRAWV(X,Y,5)
      ENDIF
C
      LINSTY=1
      CALL MVBITS(LINSTY,0,3,LINATT,10)
C
C            Draw name
C
      YL=YLEV-SIZ2+(SIZ8/2.)
      IF(FACHV.EQ.1.)YL=YLEV-SIZ8
      XL=XLEV-(SIZ4*FACHV)+SIZ8
      IF(FACHV.EQ.1.)XL=XLEV
      SZ=1.5*SIZ4*2.
      IF(FACHV.EQ.1.)SZ=1.*SIZ4*2.
      ITXCOL=5
      CALL MVBITS(ITXCOL,0,4,ITXATT,3)
      AITXCO=ITXCOL
      CALL IGSET('TXCI',AITXCO)
      CALL IGSET('TXFP',-60.)
      IWTY=IGIWTY(1)
      IF(IWTY.GT.10.OR.IWTY.LT.1)CALL IGSET('TXFP',-61.)
***      CALL IGPID(1,'Tree',NPILE,' ')
      CALL GDRAWT(XL,YL,NAME,SZ,0.,4,0)
      CALL IGSET('TXFP',2.)
      CALL IGSET('TXCI',1.)
C
      IF (IDUP.EQ.1) GO TO 20
C
C            Draw IDUP in bottom half box
C
      CALL UCTOH('<"4#>',NAMT0(1),1,5)
      IF (IDUP.LT.0) CALL UCTOH('F',NAMT0(3),1,1)
      DUP=FLOAT(IABS(IDUP))
      CALL HBCDF(DUP,8,NAMT0(6))
      CALL UCTOH('$',NAMT0(14),1,1)
      CALL UBUNCH(NAMT0,NAM0,14)
      CALL UHTOC(NAM0,4,CHNAM0,14)
      XL=XLEV+(SIZ4*FACHV)+SIZ2
      IF(FACHV.EQ.1.)XL=XLEV-(3.*SIZ8)
      YL=YLEV-SIZ4-SIZ8
      IF(FACHV.EQ.1.)YL=YLEV-SIZ2+(SIZ8/2.)
      SZ=0.7*SIZ4*2.
      IF(FACHV.EQ.1.)SZ=0.7*SIZ4*1.5
      ITXCOL=5
      CALL MVBITS(ITXCOL,0,4,ITXATT,3)
      AITXCO=ITXCOL
      CALL IGSET('TXCI',AITXCO)
      CALL IGSET('TXFP',2.)
***      CALL IGPID(1,'Tree',NPILE,' ')
      CALL GDRAWT(XL,YL,CHNAM0,SZ,0.,10,-1)
      CALL IGSET('TXCI',1.)
C
   20 CONTINUE
C
      IF (IDRAW.NE.1) GO TO 999
C
C             Draw a normal view of the volume
C
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF(IVO.LE.0) GO TO 999
C
      JVO=LQ(JVOLUM-IVO)
C
C                GET PARAMETERS
C
      CALL GFPARA(NAME,1,0,NPAR,NATT,PAR,ATT)
      IF(NPAR.LE.0) GO TO 999
C
      ISHAPE=Q(JVO+2)
C
C             CHECK parametER SIZES
C
      PARMAX=-1.
      DO 50  I=1,NPAR
         IF(ISHAPE.GT.14)GO TO 50
         IF(NNDM(ISHAPE).LE.0) GO TO 40
C
         NDM=NNDM(ISHAPE)
         DO 30  IDM=1,NDM
C
            IF(I.EQ.INDM(IDM,ISHAPE)) GO TO 50
C
   30    CONTINUE
   40    CONTINUE
C
         ABSPAR=ABS(PAR(I))
         PARMAX=MAX(PARMAX,ABSPAR)
   50 CONTINUE
      XSCALE=SIZE/(3.*PARMAX)
      YSCALE=XSCALE
      XPOS=XLEV-(SIZE*FACHV*2)+SIZ8
      YPOS=YLEV
*      YPOS=YLEV+SIZE
C
      CALL GFATT(NAME,'SEEN',KSEEN)
C
C             Add local value KSEEN 1 to starting node of tree
C
      KSEEN=KSEEN+110
      ISEEN=KSEEN
      CALL GSATT(NAME,'SEEN',ISEEN)
C
      CALL GDNSON(NAME,NSON,IDIV)
      DO 60  N=1,NSON
         CALL GDSON(N,NAME,ISON)
         CALL GFATT(ISON,'SEEN',KSEEN)
C
C             ISON is a volume with multeplicity;
C             first occurrence has already been set
C
         IF (KSEEN.GT.50) GO TO 60
C
C             Add local value SEEN -2 to each one-level-down node
C
         KSEEN=KSEEN+80
         ISEEN=KSEEN
         CALL GSATT(ISON,'SEEN',ISEEN)
   60 CONTINUE
      CALL GDRAW(NAME,45.,135.,0.,XPOS,YPOS,XSCALE,YSCALE)
C
C             Reset global SEEN values
C
      DO 70  IVO=1,NVOLUM
         CALL UHTOC(IQ(JVOLUM+IVO),4,NAMSEE,4)
         CALL GFATT(NAMSEE,'SEEN',KSEEN)
         IF (KSEEN.LT.50) GO TO 70
         ISEENL=KSEEN/10.+0.5
         ISEENG=KSEEN-ISEENL*10.
         CALL GSATT(NAMSEE,'SEEN',ISEENG)
   70 CONTINUE
C
  999 END

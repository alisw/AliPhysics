*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GETSHP(ISHAPE,PAR)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Converts GEANT SHAPES to SET. Calls  Subroutines for     *
C.    *       shape definition in SET are called
C.    *                                                                *
C.    *                                                                *
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
C.    *  13       ELTU     A,B,DZ                                      *
C.    *  14       HYPE     RMIN,RMAX,DZ,PHI
C.    *  NSGTRA   GTRA     DZ,TH,PHI,TWIST,Y1,XL1,XH1,TH1,Y2,XL2,XH2,..*
C.    *  NSCTUB   CTUB     RMIN,RMAX,DZ,PHIMIN,PHIMAX,LXL,LYL,LZL,LXH,.*
C.    *                                                                *
C.    *    ==>Called by : GDEXCA
C.    *                                                                *
C.    *       Authors : NJ Hoimyr, J Vuoskoski     *********           *
C.    ******************************************************************
C. -------------------------------------------------------
#include "geant321/gcsetf.inc"
 
 
      DIMENSION PAR(100)
C
C----------------------------------------------------------
C     IF SHAPE THEN LOOP FOR SET OUTPUT
C
      IF (ISHAPE.EQ.1) THEN
C
C       BOX
         CALL GBOSET(PAR)
C
      ELSEIF (ISHAPE.EQ.2) THEN
C
C       TRD1
         CALL GTRDSE(PAR,ISHAPE)
C
      ELSEIF (ISHAPE.EQ.3) THEN
C
C       TRD2
         CALL GTRDSE(PAR,ISHAPE)
C
      ELSEIF (ISHAPE.EQ.4) THEN
C       TRAP
         CALL GTRASE(PAR)
C
      ELSEIF (ISHAPE.EQ.5) THEN
C
C             TUBE
         CALL GTUSET(PAR,ISHAPE)
C
      ELSEIF (ISHAPE.EQ.6) THEN
C
C             TUBS
         CALL GTUSET(PAR,ISHAPE)
C
      ELSEIF (ISHAPE.EQ.7) THEN
C
C             CONE
         CALL GCOSET(PAR,ISHAPE)
C
      ELSEIF (ISHAPE.EQ.8) THEN
C
C             CONS
         CALL GCOSET(PAR,ISHAPE)
C
      ELSEIF (ISHAPE.EQ.9) THEN
C
C             SPHE
         CALL GSPSET(PAR)
C
      ELSEIF (ISHAPE.EQ.10) THEN
C
C             PARA
         CALL GPASET(PAR)
C
      ELSEIF (ISHAPE.EQ.11) THEN
C
C             PGON
         CALL GPGSET(PAR)
C
      ELSEIF (ISHAPE.EQ.12) THEN
C
C             PCON
         CALL GPCSET(PAR)
C
      ELSEIF (ISHAPE.EQ.13) THEN
C
C             ELTU
         CALL GELSET(PAR)
C
C      ELSEIF (ISHAPE.EQ.14) THEN
C     not ready yet
C
C             HYPE
C        CALL HYPSET(PAR)
C
 
      ELSEIF (ISHAPE.EQ.28) THEN
C
C             GTRA
         CALL GTRSET(PAR)
C
      ELSEIF (ISHAPE.EQ.29) THEN
C
C             CTUB
         CALL GCTSET(PAR)
      ELSE
C
C   New shapes could be included
C
      END IF
C
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:36  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*=== bdpree ===========================================================*
*
      SUBROUTINE FDPREE
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 16 september 1991 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 23-feb-93     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/nucgeo.inc"
#include "geant321/nuclev.inc"
      LOGICAL LABRS1,LELST1,LINEL1, LCHEX1, LBSRP1, LABST1
      DIMENSION JUSNU1(160,2), MAGNU1(8)
*
*
      DATA LABRS1, LELST1, LINEL1, LCHEX1, LBSRP1, LABST1 / 6*.FALSE. /
      DATA POTBA1 / POTBA0 /, POTME1 / POTME0 /, WLLRE1 / 0.D+00 /
      DATA JUSNU1 / 320 * 0 /, INUCL1 / 1 /, IEVPR1 / 0 /
      DATA MAGNU1 / 2, 8, 20, 28, 50, 82, 126, 160 /
*
      LABRST = LABRS1
      LELSTC = LELST1
      LINELS = LINEL1
      LCHEXC = LCHEX1
      LABSRP = LBSRP1
      LABSTH = LABST1
*
      POTBAR = POTBA1
      POTMES = POTME1
      WLLRES = WLLRE1
*
      DO 10 I=1,160
        JUSNUC(I,1) = JUSNU1(I,1)
        JUSNUC(I,2) = JUSNU1(I,2)
  10  CONTINUE
      INUCLV = INUCL1
      IEVPRE = IEVPR1
*
      DO 20 I=1,8
        MAGNUM(I) = MAGNU1(I)
 20   CONTINUE
*=== End of Block Data Bdpree =========================================*
      END

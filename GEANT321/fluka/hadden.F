*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE HADDEN.FOR
*COPY HADDEN
*
*=== hadden ===========================================================*
*
      SUBROUTINE HADDEN
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*  Hadden89: slight modifications by A. Ferrari, INFN-Milan and CERN-SPS
*            7-9-89
*----------------------------------------------------------------------*
*
#include "geant321/decayc2.inc"
#include "geant321/part3.inc"
#include "geant321/split.inc"
*
      COMMON /FKNAMS/ ANAME(110),ZKNAME(460),RKNAME(268)
C
      COMMON / FKABLT / AM   (110), GA   (110), TAU  (110), ICH   (110),
     &                  IBAR (110), K1   (110), K2   (110)
      CHARACTER*8 ANAME,ZKNAME,RKNAME,ANAMZ,ZKNAMZ
      COMMON /FKADHP/AMZ(16),GAZ(16),TAUZ(16),ICHZ(16),IBARZ(16),K1Z(16)
     *,K2Z(16),WTZ(153),II22,
     * NZKZ(153,3)
      COMMON/FKADHN/ANAMZ(16),ZKNAMZ(153)
      DATA IRETUR/0/
      IRETUR=IRETUR+1
      IF (IRETUR.GT.1) RETURN
      DO 1 I=1,94
         ANAME(I) = ANAMC(I)
         AM   (I) = AMC  (I)
         GA   (I) = GAC  (I)
         TAU  (I) = TAUC (I)
         ICH  (I) = ICHC (I)
         IBAR (I) = IBARC(I)
         K1   (I) = K1C  (I)
         K2   (I) = K2C  (I)
   1  CONTINUE
      AM(31) = 0.48D0
      DO 5 I=26,30
         K1(I) = 452
         K2(I) = 452
   5  CONTINUE
      DO 2 I=1,307
         ZKNAME(I) = ZKNAMC(I)
         WT    (I) = WTC   (I)
         NZK (I,1) = NZKC(I,1)
         NZK (I,2) = NZKC(I,2)
         NZK (I,3) = NZKC(I,3)
   2  CONTINUE
      DO 3 I=1,16
         L=I+94
         ANAME(L) = ANAMZ(I)
         AM   (L) = AMZ  (I)
         GA   (L) = GAZ  (I)
         TAU  (L) = TAUZ (I)
         ICH  (L) = ICHZ (I)
         IBAR (L) = IBARZ(I)
         K1   (L) = K1Z  (I)
         K2   (L) = K2Z  (I)
   3  CONTINUE
      DO 4 I=1,153
         L=I+307
         ZKNAME(L) = ZKNAMZ(I)
         WT    (L) = WTZ   (I)
         NZK (L,3) = NZKZ(I,3)
         NZK (L,2) = NZKZ(I,2)
         NZK (L,1) = NZKZ(I,1)
   4  CONTINUE
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:35  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE BDNOPT.FOR
*COPY BDNOPT
*
*=== bdnopt ===========================================================*
*==                                                                    *
      SUBROUTINE FDNOPT
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*   Created on  20 september 1989    by  Alfredo Ferrari - Infn Milan  *
*                                                                      *
*         Last change on 09-mar-93   by  Alfredo Ferrari               *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/fheavy.inc"
#include "geant321/finuc.inc"
#include "geant321/parevt.inc"
#include "geant321/resnuc.inc"
*
      CHARACTER*8 ANHEA1(12)
      DIMENSION AMHEA1(12),ICHEA1(12),IBHEA1(12)
      LOGICAL LDIFF1,LPOWE1,LINCT1,LEVPR1,LHEAV1,LDEEX1,LGDHP1
     &       ,LPREE1,LHLFI1,LPRFI1,LPARW1
      DIMENSION LDIFF1(39)
      LOGICAL LRNFS1
 
* /fheavy/
      DATA AMHEA1 / 12 * 0.D+00 /
      DATA ANHEA1 / 'NEUTRON ', 'PROTON  ', 'DEUTERON', '3-H     ',
     &              '3-He    ', '4-He    ', 'H-FRAG-1', 'H-FRAG-2',
     &              'H-FRAG-3', 'H-FRAG-4', 'H-FRAG-5', 'H-FRAG-6'/
      DATA ICHEA1 / 0, 1, 1, 1, 2, 2, 6*0 /,
     &     IBHEA1 / 1, 1, 2, 3, 3, 4, 6*0 /
      DATA NPHEA1 / 0 /
 
* /finuc/
      DATA NP1/ 0 /, TV1/ 0.D+00 /, TVCMS1/ 0.D+00 /, TVREC1 / 0.D+00/,
     &     TVHEA1 / 0.D+00 /, TVBIN1 / 0.D+00 /
 
* /parevt/
      DATA DPOWE1 /-13.D+00 /, FSPRD1 / 0.6D+00 /, FSHPF1 / 0.0D+00 /,
     &     RN1GS1 /-1.0D+00 /, RN2GS1 /-1.0D+00 /
      DATA LDIFF1 /  .TRUE., .FALSE., 6 * .TRUE., .FALSE., 8 * .TRUE.,
     &              .FALSE., 4 * .TRUE., .TRUE., 2 * .TRUE.,
     &              8 * .FALSE., .TRUE., .FALSE., .TRUE., .FALSE.,
     &              .TRUE., .FALSE. /
      DATA LPOWE1 / .TRUE.  /, LINCT1 / .TRUE.  /, LEVPR1 / .TRUE.  /,
     &     LHEAV1 / .TRUE.  /, LDEEX1 / .TRUE.  /, LGDHP1 / .TRUE.  /,
     &     LPREE1 / .TRUE.  /, LHLFI1 / .FALSE. /, LPRFI1 / .FALSE. /,
     &     LPARW1 / .TRUE. /
      DATA ILVMO1 / 0 /
 
* /resnuc/
      DATA IEVPL1 / 0 /, IEVPH1 / 0 /, IEVNE1 / 0 /, IEVPR1 / 0 /,
     &     IEVTR1 / 0 /, IEVDE1 / 0 /, IEV3H1 / 0 /, IEV4H1 / 0 /,
     &     IDEEX1 / 0 /
      DATA LRNFS1 / .FALSE. /
*
      DO 10 J=1,12
         AMHEAV(J) = AMHEA1(J)
         ANHEAV(J) = ANHEA1(J)
         ICHEAV(J) = ICHEA1(J)
         IBHEAV(J) = IBHEA1(J)
  10  CONTINUE
      NPHEAV = NPHEA1
*
      NP=NP1
      TV=TV1
      TVCMS=TVCMS1
      TVRECL=TVREC1
      TVHEAV=TVHEA1
      TVBIND=TVBIN1
*
      DPOWER = DPOWE1
      FSPRD0 = FSPRD1
      FSHPFN = FSHPF1
      RN1GSC = RN1GS1
      RN2GSC = RN2GS1
      DO 20 J=1,39
        LDIFFR(J) = LDIFF1(J)
  20  CONTINUE
      LPOWER = LPOWE1
      LINCTV = LINCT1
      LEVPRT = LEVPR1
      LHEAVY = LHEAV1
      LDEEXG = LDEEX1
      LGDHPR = LGDHP1
      LPREEX = LPREE1
      LHLFIX = LHLFI1
      LPRFIX = LPRFI1
      LPARWV = LPARW1
      ILVMOD = ILVMO1
*
      IEVAPL = IEVPL1
      IEVAPH = IEVPH1
      IEVNEU = IEVNE1
      IEVPRO = IEVPR1
      IEVTRI = IEVTR1
      IEVDEU = IEVDE1
      IEV3HE = IEV3H1
      IEV4HE = IEV4H1
      IDEEXG = IDEEX1
      LRNFSS = LRNFS1
*=                                               end*block.bdnopt      *
      END

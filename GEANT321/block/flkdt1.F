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
*$ CREATE BLKDT1.FOR
*COPY BLKDT1
*
*=== blkdt1 ===========================================================*
*==                                                                    *
      SUBROUTINE FLKDT1
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*----------------------------------------------------------------------*
*     Version May 81                 Pertti Aarnio                     *
*                                    Helsinki University of            *
*                                    Technology, Finland               *
*     Last change  17  march  92  by Alfredo Ferrari                   *
*                                                                      *
*                                                                      *
*     This is a block data subroutine of Fluka86, which should contain *
*     the common blocks, which can and should be initialised before    *
*     the run.                                                         *
*                                                                      *
*     Note that most of the default values for input parameters are    *
*     set here.                                                        *
*----------------------------------------------------------------------*
* The commons are INCLUDEd from Fluka87 Maclib                         *
*                 ---> Now from Flukadd Maclib  (A. Fasso' 1989)       *
*----------------------------------------------------------------------*
#include "geant321/paprop.inc"
*
*----------------------------------------------------------------------*
*
      CHARACTER*8 GENPA1(30), BTYPE1(39)
      LOGICAL LHADR1(39), LPRBD1, LPRBI1
      DIMENSION ISCOR1(10), IJDIS1(39)
* /Paprop/                                                             *
      DATA BTYPE1/'PROTON  ' , 'APROTON ' , 'ELECTRON' ,
     1            'POSITRON' , 'NEUTRIE ' , 'ANEUTRIE' ,
     2            'PHOTON  ' , 'NEUTRON ' , 'ANEUTRON' ,
     3            'MUON+   ' , 'MUON-   ' , 'KAONLONG' ,
     4            'PION+   ' , 'PION-   ' , 'KAON+   ' ,
     5            'KAON-   ' , 'LAMBDA  ' , 'ALAMBDA ' ,
     6            'KAONSHRT' , 'SIGMA-  ' , 'SIGMA+  ' ,
     7            'SIGMAZER' , 'PIZERO  ' , 'KAONZERO' ,
     8            'AKAONZER' , 'RESERVED' , 'BLANK   ' ,
     9            'BLANK   ' , 'BLANK   ' , 'HEAVYION' ,
     A            'ASIGMA- ' , 'ASIGMAZE' , 'ASIGMA+ ' ,
     B            'XSIZERO ' , 'AXSIZERO' , 'XSI-    ' ,
     C            'AXSI+   ' , 'OMEGA-  ' , 'AOMEGA+ ' /
      DATA ISCOR1 / 201, 9*0 /
      DATA GENPA1 / 'ALL-PART', 'ALL-CHAR', 'ALL-NEUT', 'ALL-NEGA',
     &  'ALL-POSI', 'NUCLEONS', 'NUC&PI+-', 'ENERGY  ', 'PIONS+- ',
     &  'BEAMPART', 'EM-ENRGY', 19*'BLANK   ' /
      DATA IJDIS1 / 4*0,2*1,33*0 /
      DATA LHADR1 / 2 *.TRUE., 5 * .FALSE., 2 * .TRUE., 2 * .FALSE.,
     &              28 * .TRUE. /
      DATA LPRBD1, LPRBI1 / 2*.FALSE. /
*
      DO 10 J=1,39
         BTYPE(J)=BTYPE1(J)
         IJDISC(J)=IJDIS1(J)
         LHADRO(J)=LHADR1(J)
  10  CONTINUE
      DO 20 J=1,30
         GENPAR(J)=GENPA1(J)
   20 CONTINUE
      DO 30 J=1,10
         ISCORE(J)=ISCOR1(J)
   30 CONTINUE
      LPRBSD = LPRBD1
      LPRBSI = LPRBI1
*=                                               end*block.blkdt1      *
      END

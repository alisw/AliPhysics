*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 14/09/94  14.33.42  by  S.Ravndal
*-- Author :
      SUBROUTINE GPART
C.
C.    *******************************************************************
C.    *                                                                 *
C.    *  Define standard GEANT particles plus selected decay modes      *
C.    *  and branching ratios.                                          *
C.    *                                                                 *
C.    *  All data taken from : M. AGUILAR-BENITEZ et al,                *
C.    *                        Review of Particle Properties,           *
C.    *                        Rev. Mod. Phys. 56(1984).                *
C.    *                                                                 *
C.    *  Updated to:           Phys.Rev. D50,1173(1994)                 *
C.    *                        S. Ravndal                               *
C.    *                                                                 *
C.    *                                                                 *
C.    *  In the case of W and Z, the lifetimes are calculated from      *
C.    *  quoted upper limits on widths.                                 *
C.    *                                                                 *
C.    *  NPAR           Number of parent particles defined for decay.   *
C.    *  IPAR           List of parent partilces allowed to decay.      *
C.    *                 Currently set up for pi0,pi+,pi-,K0long,K+,K-,  *
C.    *                 K0short,eta,lambda,sigma+,sigma0,sigma-,xi0,    *
C.    *                 xi-,omega-,antilambda,antisigma -,antisigma 0,  *
C.    *                 antisigma +,antixi 0,antixi +,antiomega +       *
C.    *                 decays.                                         *
C.    *  MODE(I,J)      I'th decay mode of J'th particle defined in     *
C.    *                 IPAR.                                           *
C.    *  BRATIO(I,J)    Branching ratio for I'th decay mode of J'th     *
C.    *                 particle in IPAR.                               *
C.    *                                                                 *
C.    *    ==>Called by : <USER>, UGINIT                                *
C.    *       Author    G.Patrick  *********                            *
C.    *                                                                 *
C.    *******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
C
      DIMENSION BRATIO(6,24),BRTIO1(6,20),BRTIO2(6,4)
      DIMENSION IPAR(24)
      DIMENSION MODE(6,24),MODE1(6,20),MODE2(6,4)
      DIMENSION UBUF(1)
C
      EQUIVALENCE(BRTIO1(1,1),BRATIO(1,1)),(BRTIO2(1,1),BRATIO(1,21))
      EQUIVALENCE(MODE1(1,1),MODE(1,1)),(MODE2(1,1),MODE(1,21))
C
      DATA IPAR/ 7, 8, 9,10,11,12,16,17,18,19,20,21,22,23,24,26,27,
     *          28,29,30,31,32,5,6/
      DATA BRTIO1/ 98.802, 1.198, 0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *             21.60, 19.35, 19.35, 13.50, 13.50, 12.38,
     *             63.51, 21.17,  5.59,  4.82,  3.18,  1.73,
     *             63.51, 21.17,  5.59,  4.82,  3.18,  1.73,
     *             68.61, 31.39,  0.,    0.,    0.,    0.,
     *             38.80, 31.90, 23.6,   4.88,  0.5,   0.071,
     *             63.90, 35.80,  0.,    0.,    0.,    0.,
     *             51.57, 48.30,  0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *             67.80, 23.60,  8.60,  0.,    0.,    0.,
     *             63.90, 35.80,  0.,    0.,    0.,    0.,
     *             51.57, 48.30,  0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0./
      DATA BRTIO2/100.00,  0.,    0.,    0.,    0.,    0.,
     *             67.80, 23.60,  8.60,  0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0.,
     *            100.00,  0.,    0.,    0.,    0.,    0./
      DATA MODE1/  101,30201,    0,    0,    0,    0,
     *             405,    0,    0,    0,    0,    0,
     *             406,    0,    0,    0,    0,    0,
     *             70707,40308,40209,40608,40509,70908,
     *             405,  708,90808,70402,70405,70708,
     *             406,  709,80909,70403,70406,70709,
     *             908,  707,    0,    0,    0,    0,
     *             101,70707,70908,10908,10302,10107,
     *             914,  713,    0,    0,    0,    0,
     *             714,  813,    0,    0,    0,    0,
     *             118,    0,    0,    0,    0,    0,
     *             913,    0,    0,    0,    0,    0,
     *             718,    0,    0,    0,    0,    0,
     *             918,    0,    0,    0,    0,    0,
     *            1218,  922,  723,    0,    0,    0,
     *             815,  725,    0,    0,    0,    0,
     *             715,  925,    0,    0,    0,    0,
     *             126,    0,    0,    0,    0,    0,
     *             825,    0,    0,    0,    0,    0,
     *             726,    0,    0,    0,    0,    0/
      DATA MODE2/  826,    0,    0,    0,    0,    0,
     *            1126,  830,  731,    0,    0,    0,
     *           40402,    0,    0,    0,    0,    0,
     *           40403,    0,    0,    0,    0,    0/
      DATA NPAR/24/
      DATA UBUF/0./
C.
C.    -------------------------------------------------------------------
C.
      REMASS=EMASS
      REMMU =EMMU
      RPMASS=PMASS
      CALL GSPART( 1,'GAMMA$      ',1,0.      , 0.,1.E+15,UBUF,0)
      CALL GSPART( 2,'POSITRON$   ',2,REMASS  , 1.,1.E+15,UBUF,0)
      CALL GSPART( 3,'ELECTRON$   ',2,REMASS  ,-1.,1.E+15,UBUF,0)
      CALL GSPART( 4,'NEUTRINO$   ',3,0.      , 0.,1.E+15,UBUF,0)
      CALL GSPART( 5,'MUON +$     ',5,REMMU   , 1.,2.19703E-06,UBUF,0)
      CALL GSPART( 6,'MUON -$     ',5,REMMU   ,-1.,2.19703E-06,UBUF,0)
      CALL GSPART( 7,'PION 0$     ',3,0.1349764,0.,0.84E-16,UBUF,0)
      CALL GSPART( 8,'PION +$     ',4,0.1395700, 1.,2.603E-08,UBUF,0)
      CALL GSPART( 9,'PION -$     ',4,0.1395700,-1.,2.603E-08,UBUF,0)
      CALL GSPART(10,'KAON 0 LONG$',3,0.497672, 0.,5.17E-08,UBUF,0)
      CALL GSPART(11,'KAON +$     ',4,0.493677, 1.,1.237E-08,UBUF,0)
      CALL GSPART(12,'KAON -$     ',4,0.493677,-1.,1.237E-08,UBUF,0)
      CALL GSPART(13,'NEUTRON$    ',3,0.93956563, 0.,8.870E+02,UBUF,0)
      CALL GSPART(14,'PROTON$     ',4,RPMASS  , 1.,1.000000E+15,UBUF,0)
      CALL GSPART(15,'ANTIPROTON$ ',4,RPMASS  ,-1.,1.000000E+15,UBUF,0)
      CALL GSPART(16,'KAON 0 SHORT$',3,0.497672, 0.,8.926E-11,UBUF,0)
      CALL GSPART(17,'ETA$        ',3,0.54745 , 0.,5.485E-19,UBUF,0)
      CALL GSPART(18,'LAMBDA$     ',3,1.115684, 0.,2.632E-10,UBUF,0)
      CALL GSPART(19,'SIGMA +$    ',4,1.18937 , 1.,0.799E-10,UBUF,0)
      CALL GSPART(20,'SIGMA 0$    ',3,1.19255 , 0.,7.4E-20,UBUF,0)
      CALL GSPART(21,'SIGMA -$    ',4,1.197436,-1.,1.479E-10,UBUF,0)
      CALL GSPART(22,'XI 0$       ',3,1.3149  , 0.,2.9E-10,UBUF,0)
      CALL GSPART(23,'XI -$       ',4,1.32132 ,-1.,1.639E-10,UBUF,0)
      CALL GSPART(24,'OMEGA -$    ',4,1.67245 ,-1.,0.822E-10,UBUF,0)
      CALL GSPART(25,'ANTINEUTRON$',3,0.93956563, 0.,8.870E+02,UBUF,0)
      CALL GSPART(26,'ANTILAMBDA$ ',3,1.115684, 0.,2.632E-10,UBUF,0)
      CALL GSPART(27,'ANTISIGMA -$',4,1.18937 ,-1.,0.799E-10,UBUF,0)
      CALL GSPART(28,'ANTISIGMA 0$',3,1.19255 , 0.,7.4E-20,UBUF,0)
      CALL GSPART(29,'ANTISIGMA +$',4,1.197436, 1.,1.479E-10,UBUF,0)
      CALL GSPART(30,'ANTIXI 0$   ',3,1.3149  , 0.,2.9E-10,UBUF,0)
      CALL GSPART(31,'ANTIXI +$   ',4,1.32132 , 1.,1.639E-10,UBUF,0)
      CALL GSPART(32,'ANTIOMEGA +$',4,1.67245 , 1.,0.822E-10,UBUF,0)
*     CALL GSPART(33,'TAU +$      ',4,1.7771  , 1.,3.040000E-13,UBUF,0)
*     CALL GSPART(34,'TAU -$      ',4,1.7771  ,-1.,3.040000E-13,UBUF,0)
*     CALL GSPART(35,'D +$        ',4,1.8694  , 1.,1.062000E-12,UBUF,0)
*     CALL GSPART(36,'D -$        ',4,1.8694  ,-1.,1.062000E-12,UBUF,0)
*     CALL GSPART(37,'D 0$        ',3,1.8646  , 0.,4.280000E-13,UBUF,0)
*     CALL GSPART(38,'ANTI D 0$   ',3,1.8646  , 0.,4.280000E-13,UBUF,0)
*     CALL GSPART(39,'DS+$        ',4,1.9685  , 1.,4.360000E-13,UBUF,0)
*     CALL GSPART(40,'DS-$        ',4,1.9685  ,-1.,4.360000E-13,UBUF,0)
*     CALL GSPART(41,'LAMBDA C +$ ',4,2.2851  , 1.,1.790000E-13,UBUF,0)
*     CALL GSPART(42,'W +$        ',4,80.220  , 1.,3.160000E-25,UBUF,0)
*     CALL GSPART(43,'W -$        ',4,80.220  ,-1.,2.640000E-25,UBUF,0)
*     CALL GSPART(44,'Z 0$        ',3,91.187  , 0.,7.740000E-26,UBUF,0)
      CALL GSPART(45,'DEUTERON$   ',8,1.875613,+1.,1.000000E+15,UBUF,0)
      CALL GSPART(46,'TRITON$     ',8,2.80925 ,+1.,1.000000E+15,UBUF,0)
      CALL GSPART(47,'ALPHA$      ',8,3.727417,+2.,1.000000E+15,UBUF,0)
      CALL GSPART(48,'GEANTINO$   ',6,0.      , 0.,1.000000E+15,UBUF,0)
      CALL GSPART(49,'HE3$        ',8,2.80923 ,+2.,1.000000E+15,UBUF,0)
      CALL GSPART(50,'Cerenkov$   ',7,0.      , 0.,1.000000E+15,UBUF,0)
C
C     Define decay modes.
C
      DO 10 I=1,NPAR
      CALL GSDK(IPAR(I),BRATIO(1,I),MODE(1,I))
   10 CONTINUE
C
  99  RETURN
      END

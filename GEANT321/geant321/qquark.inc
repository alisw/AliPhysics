*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_QQUARK_INC
#define CERNLIB_GEANT321_QQUARK_INC
*
*
* qquark.inc
*
*
*=== qquark ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on    6 february 1991    by        Alfredo Ferrari       *
*                                                  INFN - Milan        *
*                                                                      *
*     Last change  on  6 february 1991  by       Alfredo Ferrari       *
*                                                                      *
*     Included in the following routines :                             *
*                                                                      *
*                     COREVT                                           *
*                     CORRIN                                           *
*                     HADEVV                                           *
*                     HADEVT                                           *
*                     NUCEVV                                           *
*                     NUCEVT                                           *
*                                                                      *
*     Quark content of particles:                                      *
*          index   quark   el. charge  bar. charge  isospin  isospin3  *
*              1 = u          2/3          1/3        1/2       1/2    *
*             -1 = ubar      -2/3         -1/3        1/2      -1/2    *
*              2 = d         -1/3          1/3        1/2      -1/2    *
*             -2 = dbar       1/3         -1/3        1/2       1/2    *
*              3 = s         -1/3          1/3         0         0     *
*             -3 = sbar       1/3         -1/3         0         0     *
*              4 = c          2/3          1/3         0         0     *
*             -4 = cbar      -2/3         -1/3         0         0     *
*              5 = b         -1/3          1/3         0         0     *
*             -5 = bbar       1/3         -1/3         0         0     *
*              6 = t          2/3          1/3         0         0     *
*             -6 = tbar      -2/3         -1/3         0         0     *
*                                                                      *
*         Mquark = particle quark composition (Paprop numbering)       *
*         Iqechr = electric charge ( in 1/3 unit )                     *
*         Iqbchr = baryonic charge ( in 1/3 unit )                     *
*         Iqichr = isospin ( in 1/2 unit ), z component                *
*         Iqschr = strangeness                                         *
*         Iqcchr = charm                                               *
*         Iquchr = beauty                                              *
*         Iqtchr = ......                                              *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON / FKQUAR / IQECHR (-6:6), IQBCHR (-6:6), IQICHR (-6:6),
     &                  IQSCHR (-6:6), IQCCHR (-6:6), IQUCHR (-6:6),
     &                  IQTCHR (-6:6), MQUARK (3,39)
 

#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:35  cernlib
* Version 7.04
*
* Revision 1.2  1996/10/28 15:02:12  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       SUBROUTINE GRVGAH0 (ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*      G R V - P H O T O N - P A R A M E T R I Z A T I O N S      *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE :                *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/31             *
*                                                                 *
*    THE OUTPUT IS ALWAYS   1./ ALPHA(EM) * X * PARTON DENSITY    *
*    output modified by HPB to be always    X * PARTON DENSITY    *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   *
*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
*   / HO) AND  1.E6 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
*                                                                 *
*              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              *
*         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         *
*                                                                 *
*      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     *
*      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     *
*                                                                 *
*      HO DISTRIBUTIONS REFER TO THE DIS(GAMMA) SCHEME, SEE :     *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/26             *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
       IMPLICIT REAL (A - Y)
#include "pdf/expdp.inc"
     +        ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL
       REAL  X, Q
       DATA ALPHEM/7.29927D-3/
       X = ZX
       Q = ZQ
       MU2  = 0.3
       LAM2 = 0.248 * 0.248
       Q2 = Q*Q
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  1.447
       BE =  0.848
       AK =  0.527 + 0.200 * S  - 0.107 * S2
       BK =  7.106 - 0.310 * SS - 0.786 * S2
       AG =  0.197 + 0.533 * S
       BG =  0.062 - 0.398 * S  + 0.109 * S2
       C  =          0.755 * S  - 0.112 * S2
       D  =  0.318 - 0.059 * S
       E  =  4.225 + 1.708 * S
       ES =  1.752 + 0.866 * S
       U0 =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZUV = U0 * ALPHEM
       ZUB = ZUV
C...X * D = X * DBAR :
       AL =  1.424
       BE =  0.770
       AK =  0.500 + 0.067 * SS - 0.055 * S2
       BK =  0.376 - 0.453 * SS + 0.405 * S2
       AG =  0.156 + 0.184 * S
       BG =   0.0  - 0.528 * S  + 0.146 * S2
       C  =  0.121 + 0.092 * S
       D  =  0.379 - 0.301 * S  + 0.081 * S2
       E  =  4.346 + 1.638 * S
       ES =  1.645 + 1.016 * S
       D0  =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZDV = D0 * ALPHEM
       ZDB = ZDV
C...X * G :
       AL =  0.661
       BE =  0.793
       AK =  0.537 - 0.600 * SS
       BK =  6.389              - 0.953 * S2
       AG =  0.558 - 0.383 * SS + 0.261 * S2
       BG =   0.0  - 0.305 * S
       C  = -0.222              + 0.078 * S2
       D  =  0.153 + 0.978 * S  - 0.209 * S2
       E  =  1.429 + 1.772 * S
       ES =  3.331 + 0.806 * S
       G0 =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZGL = G0 * ALPHEM
C...X * S = X * SBAR :
       SF =   0.0
       AL =  1.578
       BE =  0.863
       AK =  0.622 + 0.332 * S  - 0.300 * S2
       BK =  2.469
       AG =  0.211 - 0.064 * SS - 0.018 * S2
       BG = -0.215 + 0.122 * S
       C  =  0.153
       D  =   0.0  + 0.253 * S  - 0.081 * S2
       E  =  3.990 + 2.014 * S
       ES =  1.720 + 0.986 * S
       S0 =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZSB = S0 * ALPHEM
C...X * C = X * CBAR :
       SF =  0.820
       AL =  0.929
       BE =  0.381
       AK =  1.228 - 0.231 * S
       BK =  3.806             - 0.337 * S2
       AG =  0.932 + 0.150 * S
       BG = -0.906
       C  =  1.133
       D  =   0.0  + 0.138 * S  - 0.028 * S2
       E  =  5.588 + 0.628 * S
       ES =  2.665 + 1.054 * S
       C0 =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZCB = C0 * ALPHEM
C...X * B = X * BBAR :
       SF =  1.297
       AL =  0.970
       BE =  0.207
       AK =  1.719 - 0.292 * S
       BK =  0.928 + 0.096 * S
       AG =  0.845 + 0.178 * S
       BG = -2.310
       C  =  1.558
       D  = -0.191 + 0.151 * S
       E  =  6.089 + 0.282 * S
       ES =  3.379 + 1.062 * S
       B0 =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZBB = B0 * ALPHEM
C
       RETURN
       END

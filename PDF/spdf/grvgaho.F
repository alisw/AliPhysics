*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:40  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       SUBROUTINE GRVGAHO (ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL)
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
       DATA ALPHEM/7.29927D-3/
       REAL  X, Q
       X = ZX
       Q = ZQ
       MU2  = 0.3
       LAM2 = 0.248 * 0.248
       Q2 = Q*Q
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  0.583
       BE =  0.688
       AK =  0.449 - 0.025 * S  - 0.071 * S2
       BK =  5.060 - 1.116 * SS
       AG =  0.103
       BG =  0.319 + 0.422 * S
       C  =  1.508 + 4.792 * S  - 1.963 * S2
       D  =  1.075 + 0.222 * SS - 0.193 * S2
       E  =  4.147 + 1.131 * S
       ES =  1.661 + 0.874 * S
       UH =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZUV = UH * ALPHEM
       ZUB = ZUV
C...X * D = X * DBAR :
       AL =  0.591
       BE =  0.698
       AK =  0.442 - 0.132 * S  - 0.058 * S2
       BK =  5.437 - 1.916 * SS
       AG =  0.099
       BG =  0.311 - 0.059 * S
       C  =  0.800 + 0.078 * S  - 0.100 * S2
       D  =  0.862 + 0.294 * SS - 0.184 * S2
       E  =  4.202 + 1.352 * S
       ES =  1.841 + 0.990 * S
       DH  =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZDV = DH * ALPHEM
       ZDB = ZDV
C...X * G :
       AL =  1.161
       BE =  1.591
       AK =  0.530 - 0.742 * SS + 0.025 * S2
       BK =  5.662
       AG =  0.533 - 0.281 * SS + 0.218 * S2
       BG =  0.025 - 0.518 * S  + 0.156 * S2
       C  = -0.282              + 0.209 * S2
       D  =  0.107 + 1.058 * S  - 0.218 * S2
       E  =   0.0  + 2.704 * S
       ES =  3.071 - 0.378 * S
       GH =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZGL = GH * ALPHEM
C...X * S = X * SBAR :
       SF =   0.0
       AL =  0.635
       BE =  0.456
       AK =  1.770 - 0.735 * SS - 0.079 * S2
       BK =  3.832
       AG =  0.084 - 0.023 * S
       BG =  0.136
       C  =  2.119 - 0.942 * S  + 0.063 * S2
       D  =  1.271 + 0.076 * S  - 0.190 * S2
       E  =  4.604 + 0.737 * S
       ES =  1.641 + 0.976 * S
       SH =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZSB = SH * ALPHEM
C...X * C = X * CBAR :
       SF =  0.820
       AL =  0.926
       BE =  0.152
       AK =  1.142 - 0.175 * S
       BK =  3.276
       AG =  0.504 + 0.317 * S
       BG = -0.433
       C  =  3.334
       D  =  0.398 + 0.326 * S  - 0.107 * S2
       E  =  5.493 + 0.408 * S
       ES =  2.426 + 1.277 * S
       CH =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZCB = CH * ALPHEM
C...X * B = X * BBAR :
       SF =  1.297
       AL =  0.969
       BE =  0.266
       AK =  1.953 - 0.391 * S
       BK =  1.657 - 0.161 * S
       AG =  1.076 + 0.034 * S
       BG = -2.015
       C  =  1.662
       D  =  0.353 + 0.016 * S
       E  =  5.713 + 0.249 * S
       ES =  3.456 + 0.673 * S
       BH =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZBB = BH * ALPHEM
c
       RETURN
       END

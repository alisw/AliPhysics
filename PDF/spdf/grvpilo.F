*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:44  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
       SUBROUTINE GRVPILO (ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*         G R V - P I O N - P A R A M E T R I Z A T I O N S       *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE :                *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/16             *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   *
*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
*   / HO) AND  1.E8 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
*   REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION IS NEG-   *
*   LIGIBLE, I.E. BELOW ABOUT 1.E-4, WERE EXCLUDED FROM THE FIT.  *
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
*   HO DISTRIBUTION REFER TO THE MS-BAR SCHEME OF BARDEEN ET AL.  *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
       IMPLICIT REAL (A - Y)
#include "pdf/expdp.inc"
     +        ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL
       REAL  X, Q
       X = ZX
       Q = ZQ
       MU2  = 0.25
       LAM2 = 0.232 * 0.232
       Q2 = Q*Q
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
C...X * VALENCE :
       NV  =  0.519 + 0.180 * S - 0.011 * S2
       AKV =  0.499 - 0.027 * S
       AGV =  0.381 - 0.419 * S
       DV  =  0.367 + 0.563 * S
       VAP =  GRVFVP (X, NV, AKV, AGV, DV)
       ZUV = VAP
       ZDV = ZUV
C...X * GLUON :
       ALG =  0.599
       BEG =  1.263
       AKG =  0.482 + 0.341 * DS
       BKG =   0.0
       AGG =  0.678 + 0.877 * S  - 0.175 * S2
       BGG =  0.338 - 1.597 * S
       CG  =   0.0  - 0.233 * S  + 0.406 * S2
       DG  =  0.390 + 1.053 * S
       EG  =  0.618 + 2.070 * S
       ESG =  3.676
       GLP =  GRVFGP(X,S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
       ZGL = GLP
C...X * QBAR (SU(3)-SYMMETRIC SEA) :
       SL  =   0.0
       ALS =   0.55
       BES =   0.56
       AKS =  2.538 - 0.763 * S
       AGS = -0.748
       BS  =  0.313 + 0.935 * S
       DS  =  3.359
       EST =  4.433 + 1.301 * S
       ESS =   9.30 - 0.887 * S
       QBP =  GRVFQBP (X, S, SL, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
       ZUDB = QBP
       ZSB = ZUDB
C...X * CBAR = X * C :
       SC  =  0.888
       ALC =   1.02
       BEC =   0.39
       AKC =   0.0
       AGC =   0.0
       BC  =  1.008
       DC  =  1.208 + 0.771 * S
       EC  =   4.40 + 1.493 * S
       ESC =  2.032 + 1.901 * S
       CBP =  GRVFQBP (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
       ZCB = CBP
C...X * BBAR = X * B :
       SBO =  1.351
       ALB =   1.03
       BEB =   0.39
       AKB =   0.0
       AGB =   0.0
       BBO =   0.0
       DB  =  0.697 + 0.855 * S
       EB  =   4.51 + 1.490 * S
       ESB =  3.056 + 1.694 * S
       BBP =  GRVFQBP (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
       ZBB = BBP
C...X * TBAR = X * T :
       TBP = 0.
       ZTB = TBP
       RETURN
       END

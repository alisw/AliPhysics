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
       SUBROUTINE GRVPIHO (ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL)
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
       MU2  = 0.3
       LAM2 = 0.248 * 0.248
       Q2 = Q*Q
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
C...X * VALENCE :
       NV  =  0.456 + 0.150 * DS + 0.112 * S - 0.019 * S2
       AKV =  0.505 - 0.033 * S
       AGV =  0.748 - 0.669 * DS - 0.133 * S
       DV  =  0.365 + 0.197 * DS + 0.394 * S
       VAP =  GRVFVP (X, NV, AKV, AGV, DV)
       ZUV = VAP
       ZDV = ZUV
C...X * GLUON :
       ALG =  1.096
       BEG =  1.371
       AKG =  0.437 - 0.689 * DS
       BKG = -0.631
       AGG =  1.324 - 0.441 * DS - 0.130 * S
       BGG = -0.955 + 0.259 * S
       CG  =  1.075 - 0.302 * S
       DG  =  1.158 + 1.229 * S
       EG  =   0.0  + 2.510 * S
       ESG =  2.604 + 0.165 * S
       GLP =  GRVFGP(X,S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
       ZGL = GLP
C...X * QBAR (SU(3)-SYMMETRIC SEA) :
       SL  =   0.0
       ALS =   0.85
       BES =   0.96
       AKS = -0.350 + 0.806 * S
       AGS = -1.663
       BS  =  3.148
       DS  =  2.273 + 1.438 * S
       EST =  3.214 + 1.545 * S
       ESS =  1.341 + 1.938 * S
       QBP =  GRVFQBP (X, S, SL, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
       ZUDB = QBP
       ZSB = ZUDB
C...X * CBAR = X * C :
       SC  =  0.820
       ALC =   0.98
       BEC =   0.0
       AKC =   0.0  - 0.457 * S
       AGC =   0.0
       BC  =  -1.00 +  1.40 * S
       DC  =  1.318 + 0.584 * S
       EC  =   4.45 + 1.235 * S
       ESC =  1.496 + 1.010 * S
       CBP =  GRVFQBP (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
       ZCB = CBP
C...X * BBAR = X * B :
       SBO =  1.297
       ALB =   0.99
       BEB =   0.0
       AKB =   0.0  - 0.172 * S
       AGB =   0.0
       BBO =   0.0
       DB  =  1.447 + 0.485 * S
       EB  =   4.79 + 1.164 * S
       ESB =  1.724 + 2.121 * S
       BBP =  GRVFQBP (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
       ZBB = BBP
C...X * TBAR = X * T :
       TBP = 0.
       ZTB = TBP
       RETURN
       END

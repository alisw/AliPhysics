      subroutine GRVP1evolve(xin,qin,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,pdf(-6:6),xval(45),qcdl4,qcdl5      
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      
      save 

      call grvpiho(xin,Qin,upv,dnv,usea,str,chm,bot,top,glu)
            
      pdf(-6)= top
      pdf(6)= top
      pdf(-5)= bot
      pdf(5 )= bot
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv+usea
      pdf(-1)= usea
      pdf(1 )= dnv+usea
      pdf(0 )= glu

      return
c      
      entry GRVP0evolve(xin,qin,pdf)

      call grvpilo(xin,Qin,upv,dnv,usea,str,chm,bot,top,glu)
            
      pdf(-6)= top
      pdf(6)= top
      pdf(-5)= bot
      pdf(5 )= bot
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv+usea
      pdf(-1)= usea
      pdf(1 )= dnv+usea
      pdf(0 )= glu

      return
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVPread(nset)
      read(1,*)nmem(nset),ndef(nset)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVPalfa(alfas,qalfa)
        call getnset(iset)
	call getnmem(iset,imem)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVPinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVPpdf(mem)
c      imem = mem
      call getnset(iset)
      call setnmem(iset,mem)
      return
c
 1000 format(5e13.5)
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:28:44  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
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
       double precision
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
c=================================================================
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:28:44  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
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
       double precision
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
c====================================================================
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:28:38  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
       FUNCTION GRVFVP (X, N, AK, AG, D)
       IMPLICIT REAL (A - Z)
       DX = SQRT (X)
       GRVFVP = N * X**AK * (1.+ AG*DX) * (1.- X)**D
       RETURN
       END
c====================================================================
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:28:37  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
       FUNCTION GRVFQBP (X, S, ST, AL, BE, AK, AG, B, D, E, ES)
       IMPLICIT REAL (A - Z)
       DX = SQRT (X)
       LX = ALOG (1./X)
       IF (S .LE. ST) THEN
          GRVFQBP = 0.0
       ELSE
          GRVFQBP = (S-ST)**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1           * EXP (-E + SQRT (ES * S**BE * LX))
       END IF
       RETURN
       END
c====================================================================
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:28:36  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
       FUNCTION GRVFGP (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT REAL (A - Z)
       DX = SQRT (X)
       LX = ALOG (1./X)
       GRVFGP = (X**AK * (AG + BG*DX + C*X) * LX**BK + S**AL
     1       * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
c====================================================================

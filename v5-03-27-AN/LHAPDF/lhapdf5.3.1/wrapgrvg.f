      subroutine GRVGevolvep0(xin,qin,p2in,ip2in,pdf)
      include 'parmsetup.inc'
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu,zbot,zchm
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      
      save 
      
      call getnset(iset)
      call getnmem(iset,imem)
      if(imem.eq.1) then
        call GRVGALO (xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)
      elseif(imem.eq.2.or.imem.eq.0) then
        q2in = qin*qin
c calls GRVGALO for charm and bottom, rest from GRSGALO
        call GRVGALO(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)
        call GRSGALO(xin,q2in,p2in,
     +               upv,dnv,usea,dsea,str,zchm,zbot,glu)
      else
        CONTINUE
      endif     
      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= bot
      pdf(5 )= bot
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv
      pdf(-1)= dsea
      pdf(1 )= dnv
      pdf(0 )= glu
      
      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVGevolvep1(xin,qin,p2in,ip2in,pdf)

      if(imem.eq.1) then
        call GRVGAH0 (xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)
      elseif(imem.eq.2 .or. imem.eq.0) then
        call GRVGAHO (xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu)
      else
        CONTINUE
      endif     

      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= bot
      pdf(5 )= bot
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv
      pdf(-1)= dsea
      pdf(1 )= dnv
      pdf(0 )= glu
      
      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVGread(nset)
      read(1,*)nmem(nset),ndef(nset)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVGalfa(alfas,qalfa)
        call getnset(iset)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVGinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVGpdf(mem)
      call getnset(iset)
      call setnmem(iset,mem)

c      imem = mem
      return
c
 1000 format(5e13.5)
      end
c
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
       double precision
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
       double precision
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE GRVGALO (ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL)
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
       double precision
     +        ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZGL
       REAL  X, Q
       DATA ALPHEM/7.29927D-3/
       X = ZX
       Q = ZQ
       MU2  = 0.25
       LAM2 = 0.232 * 0.232
       Q2 = Q*Q
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  1.717
       BE =  0.641
       AK =  0.500 - 0.176 * S
       BK = 15.00  - 5.687 * SS - 0.552 * S2
       AG =  0.235 + 0.046 * SS
       BG =  0.082 - 0.051 * S  + 0.168 * S2
       C  =   0.0  + 0.459 * S
       D  =  0.354 - 0.061 * S
       E  =  4.899 + 1.678 * S
       ES =  2.046 + 1.389 * S
       UL =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZUV = UL * ALPHEM
       ZUB = ZUV
C...X * D = X * DBAR :
       AL =  1.549
       BE =  0.782
       AK =  0.496 + 0.026 * S
       BK =  0.685 - 0.580 * SS + 0.608 * S2
       AG =  0.233 + 0.302 * S
       BG =   0.0  - 0.818 * S  + 0.198 * S2
       C  =  0.114 + 0.154 * S
       D  =  0.405 - 0.195 * S  + 0.046 * S2
       E  =  4.807 + 1.226 * S
       ES =  2.166 + 0.664 * S
       DL  =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZDV = DL * ALPHEM
       ZDB = ZDV
C...X * G :
       AL =  0.676
       BE =  1.089
       AK =  0.462 - 0.524 * SS
       BK =  5.451              - 0.804 * S2
       AG =  0.535 - 0.504 * SS + 0.288 * S2
       BG =  0.364 - 0.520 * S
       C  = -0.323              + 0.115 * S2
       D  =  0.233 + 0.790 * S  - 0.139 * S2
       E  =  0.893 + 1.968 * S
       ES =  3.432 + 0.392 * S
       GL =  GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZGL = GL * ALPHEM
C...X * S = X * SBAR :
       SF =   0.0
       AL =  1.609
       BE =  0.962
       AK =  0.470              - 0.099 * S2
       BK =  3.246
       AG =  0.121 - 0.068 * SS
       BG = -0.090 + 0.074 * S
       C  =  0.062 + 0.034 * S
       D  =   0.0  + 0.226 * S  - 0.060 * S2
       E  =  4.288 + 1.707 * S
       ES =  2.122 + 0.656 * S
       SL =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZSB = SL * ALPHEM
C...X * C = X * CBAR :
       SF =  0.888
       AL =  0.970
       BE =  0.545
       AK =  1.254 - 0.251 * S
       BK =  3.932              - 0.327 * S2
       AG =  0.658 + 0.202 * S
       BG = -0.699
       C  =  0.965
       D  =   0.0  + 0.141 * S  - 0.027 * S2
       E  =  4.911 + 0.969 * S
       ES =  2.796 + 0.952 * S
       CL =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZCB = CL * ALPHEM
C...X * B = X * BBAR :
       SF =  1.351
       AL =  1.016
       BE =  0.338
       AK =  1.961 - 0.370 * S
       BK =  0.923 + 0.119 * S
       AG =  0.815 + 0.207 * S
       BG = -2.275
       C  =  1.480
       D  = -0.223 + 0.173 * S
       E  =  5.426 + 0.623 * S
       ES =  3.819 + 0.901 * S
       BL =  GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       ZBB = BL * ALPHEM
C
       RETURN
       END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C-----------------------------------------------------------------------
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*           G R S - LO - VIRTUAL PHOTON PARAMETRIZATIONS          *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                M. GLUECK, E.REYA, M. STRATMANN :                *
*                    PHYS. REV. D51 (1995) 3220                   *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.6   AND  5.E4                 *
*                       AND (!)  Q**2 > 5 P**2                    *
*        P**2 / GEV**2  BETWEEN   0.0   AND  10.                  *
*                       P**2 = 0  <=> REAL PHOTON                 *
*             X         BETWEEN  1.E-4  AND   1.                  *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*                                                                 *
*   PLEASE REPORT ANY STRANGE BEHAVIOUR TO :                      *
*          STRAT@HAL1.PHYSIK.UNI-DORTMUND.DE                      *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS : 
*
*    X   = MOMENTUM FRACTION 
*    Q2  = SCALE Q**2 IN GEV**2
*    P2  = VIRTUALITY OF THE PHOTON IN GEV**2 
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION DIVIDED BY ALPHA_EM) :
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION) :  modified H.P.-B. 10.9.1996
*
********************************************************
      SUBROUTINE GRSGALO(DX,DQ2,DP2,
     +                     DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL)
C      subroutine grsgalo(x,q2,p2,ugam,dgam,sgam,ggam)
      implicit real*8 (a-h,o-z)
      double precision
     +          x, q2, p2, mu2, lam2,
     +          ugam, dgam, sgam, ggam,
     +          DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DGL
C
      dimension u1(40),ds1(40),g1(40)
      dimension ud2(20),s2(20),g2(20)
      dimension up0(20),dsp0(20),gp0(20)
      DATA ALPHEM/7.29927D-3/ 
c
      data u1/-0.139d0,0.783d0,0.132d0,0.087d0,0.003d0,-0.0134d0,
     +   0.009d0,-0.017d0,0.092d0,-0.516d0,-0.085d0,0.439d0,
     +   0.013d0,0.108d0,-0.019d0,-0.272d0,-0.167d0,0.138d0,
     +   0.076d0,0.026d0,-0.013d0,0.27d0,0.107d0,-0.097d0,0.04d0,
     +   0.064d0,0.011d0,0.002d0,0.057d0,-0.057d0,0.162d0,
     +   -0.172d0,0.124d0,-0.016d0,-0.065d0,0.044d0,-1.009d0,
     +   0.622d0,0.227d0,-0.184d0/
      data ds1/0.033d0,0.007d0,-0.0516d0,0.12d0,0.001d0,-0.013d0,
     +   0.018d0,-0.028d0,0.102d0,-0.595d0,-0.114d0,0.669d0,
     +   0.022d0,0.001d0,-0.003d0,-0.0583d0,-0.041d0,0.035d0,
     +   0.009d0,0.009d0,0.004d0,0.054d0,0.025d0,-0.02d0,
     +   0.007d0,0.021d0,0.01d0,0.004d0,-0.067d0,0.06d0,-0.148d0,
     +   0.13d0,0.032d0,-0.009d0,-0.06d0,0.036d0,-0.39d0,0.033d0,
     +   0.245d0,-0.171d0/
      data g1/0.025d0,0.d0,-0.018d0,0.112d0,-0.025d0,0.177d0, 
     +   -0.022d0,0.024d0,0.001d0,-0.0104d0,0.d0,0.d0,-1.082d0,
     +   -1.666d0,0.d0,0.086d0,0.d0,0.053d0,0.005d0,-0.058d0,
     +   0.034d0,0.073d0,1.08d0,1.63d0,-0.0256d0,-0.088d0,0.d0,
     +   0.d0,-0.004d0,0.016d0,0.007d0,-0.012d0,0.01d0,-0.673d0,
     +   0.126d0,-0.167d0,0.032d0,-0.227d0,0.086d0,-0.159d0/
      data ud2/0.756d0,0.187d0,0.109d0,-0.163d0,0.002d0,0.004d0,
     +   0.054d0,-0.039d0,22.53d0,-21.02d0,5.608d0,0.332d0,
     +   -0.008d0,-0.021d0,0.381d0,0.572d0,4.774d0,1.436d0,
     +   -0.614d0,3.548d0/ 
      data s2/0.902d0,0.182d0,0.271d0,-0.346d0,0.017d0,-0.01d0,
     +   -0.011d0,0.0065d0,17.1d0,-13.29d0,6.519d0,0.031d0,
     +   -0.0176d0,0.003d0,1.243d0,0.804d0,4.709d0,1.499d0,
     +   -0.48d0,3.401d0/
      data g2/0.364d0,1.31d0,0.86d0,-0.254d0,0.611d0,0.008d0,
     +   -0.097d0,-2.412d0,-0.843d0,2.248d0,-0.201d0,1.33d0,
     +   0.572d0,0.44d0,1.233d0,0.009d0,0.954d0,1.862d0,3.791d0,
     +   -0.079d0/
      data up0/1.551d0,0.105d0,1.089d0,-0.172d0,3.822d0,-2.162d0,
     +   0.533d0,-0.467d0,-0.412d0,0.2d0,0.377d0,0.299d0,0.487d0,
     +   0.0766d0,0.119d0,0.063d0,7.605d0,0.234d0,-0.567d0,
     +   2.294d0/ 
      data dsp0/2.484d0,1.214d0,1.088d0,-0.1735d0,4.293d0,
     +   -2.802d0,0.5975d0,-0.1193d0,-0.0872d0,0.0418d0,0.128d0,
     +   0.0337d0,0.127d0,0.0135d0,0.14d0,0.0423d0,6.946d0,
     +   0.814d0,1.531d0,0.124d0/
      data gp0/1.682d0,1.1d0,0.5888d0,-0.4714d0,0.5362d0,0.0127d0,
     +   -2.438d0,0.03399d0,0.07825d0,0.05842d0,0.08393d0,2.348d0,
     +   -0.07182d0,1.084d0,0.3098d0,-0.07514d0,3.327d0,1.1d0,
     +   2.264d0,0.2675d0/
c
      save u1,ds1,g1,ud2,s2,g2,up0,dsp0,gp0
c
      x  = DX
      q  = SQRT(DQ2)
      q2 = DQ2
      p2 = DP2
      mu2=0.25d0
      lam2=0.232d0*0.232d0
c
      if(p2.le.0.25d0) then
         s=log(log(q2/lam2)/log(mu2/lam2))
         lp1=0.d0
         lp2=0.d0
      else
         if(q2.lt.p2) then
            write(*,1000)
 1000       format
     +      ('  WARNING: GRSGALO has been called with Q2 < P2 !',/,
     +       '           GRSGALO is about to blow up, therefore',/,
     +       '           Q2 is set equal to P2')
            q2=p2
         endif
         s=log(log(q2/lam2)/log(p2/lam2))
         lp1=log(p2/mu2)*log(p2/mu2)
         lp2=log(p2/mu2+log(p2/mu2))
      endif
c
      alp=up0(1)+lp1*u1(1)+lp2*u1(2)
      bet=up0(2)+lp1*u1(3)+lp2*u1(4)
      a=up0(3)+lp1*u1(5)+lp2*u1(6)+
     +  (up0(4)+lp1*u1(7)+lp2*u1(8))*s
      b=up0(5)+lp1*u1(9)+lp2*u1(10)+
     +  (up0(6)+lp1*u1(11)+lp2*u1(12))*s**0.5+
     +  (up0(7)+lp1*u1(13)+lp2*u1(14))*s**2
      gb=up0(8)+lp1*u1(15)+lp2*u1(16)+
     +  (up0(9)+lp1*u1(17)+lp2*u1(18))*s+
     +  (up0(10)+lp1*u1(19)+lp2*u1(20))*s**2
      ga=up0(11)+lp1*u1(21)+lp2*u1(22)+
     +  (up0(12)+lp1*u1(23)+lp2*u1(24))*s**0.5
      gc=up0(13)+lp1*u1(25)+lp2*u1(33)+
     +  (up0(14)+lp1*u1(26)+lp2*u1(34))*s
      gd=up0(15)+lp1*u1(27)+lp2*u1(35)+
     +  (up0(16)+lp1*u1(28)+lp2*u1(36))*s
      ge=up0(17)+lp1*u1(29)+lp2*u1(37)+
     +  (up0(18)+lp1*u1(30)+lp2*u1(38))*s
      gep=up0(19)+lp1*u1(31)+lp2*u1(39)+
     +  (up0(20)+lp1*u1(32)+lp2*u1(40))*s
      upart1=grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
c
      alp=dsp0(1)+lp1*ds1(1)+lp2*ds1(2)
      bet=dsp0(2)+lp1*ds1(3)+lp2*ds1(4)
      a=dsp0(3)+lp1*ds1(5)+lp2*ds1(6)+
     +  (dsp0(4)+lp1*ds1(7)+lp2*ds1(8))*s
      b=dsp0(5)+lp1*ds1(9)+lp2*ds1(10)+
     +  (dsp0(6)+lp1*ds1(11)+lp2*ds1(12))*s**0.5+
     +  (dsp0(7)+lp1*ds1(13)+lp2*ds1(14))*s**2
      gb=dsp0(8)+lp1*ds1(15)+lp2*ds1(16)+
     +  (dsp0(9)+lp1*ds1(17)+lp2*ds1(18))*s+
     +  (dsp0(10)+lp1*ds1(19)+lp2*ds1(20))*s**2
      ga=dsp0(11)+lp1*ds1(21)+lp2*ds1(22)+
     +  (dsp0(12)+lp1*ds1(23)+lp2*ds1(24))*s
      gc=dsp0(13)+lp1*ds1(25)+lp2*ds1(33)+
     +  (dsp0(14)+lp1*ds1(26)+lp2*ds1(34))*s
      gd=dsp0(15)+lp1*ds1(27)+lp2*ds1(35)+
     +  (dsp0(16)+lp1*ds1(28)+lp2*ds1(36))*s
      ge=dsp0(17)+lp1*ds1(29)+lp2*ds1(37)+
     +  (dsp0(18)+lp1*ds1(30)+lp2*ds1(38))*s
      gep=dsp0(19)+lp1*ds1(31)+lp2*ds1(39)+
     +  (dsp0(20)+lp1*ds1(32)+lp2*ds1(40))*s
      dspart1=grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
c
      alp=gp0(1)+lp1*g1(1)+lp2*g1(2)
      bet=gp0(2)+lp1*g1(3)+lp2*g1(4)
      a=gp0(3)+lp1*g1(5)+lp2*g1(6)+
     +  (gp0(4)+lp1*g1(7)+lp2*g1(8))*s**0.5
      b=gp0(5)+lp1*g1(9)+lp2*g1(10)+
     +  (gp0(6)+lp1*g1(11)+lp2*g1(12))*s**2
      gb=gp0(7)+lp1*g1(13)+lp2*g1(14)+
     +  (gp0(8)+lp1*g1(15)+lp2*g1(16))*s
      ga=gp0(9)+lp1*g1(17)+lp2*g1(18)+
     +  (gp0(10)+lp1*g1(19)+lp2*g1(20))*s**0.5+
     +  (gp0(11)+lp1*g1(21)+lp2*g1(22))*s**2
      gc=gp0(12)+lp1*g1(23)+lp2*g1(24)+
     +  (gp0(13)+lp1*g1(25)+lp2*g1(26))*s**2
      gd=gp0(14)+lp1*g1(27)+lp2*g1(28)+
     +  (gp0(15)+lp1*g1(29)+lp2*g1(30))*s+
     +  (gp0(16)+lp1*g1(31)+lp2*g1(32))*s**2
      ge=gp0(17)+lp1*g1(33)+lp2*g1(34)+
     +  (gp0(18)+lp1*g1(35)+lp2*g1(36))*s
      gep=gp0(19)+lp1*g1(37)+lp2*g1(38)+
     +  (gp0(20)+lp1*g1(39)+lp2*g1(40))*s
      gpart1=grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
c
      s=log(log(q2/lam2)/log(mu2/lam2))
      suppr=1.d0/(1.d0+p2/0.59d0)**2
c
      alp=ud2(1)
      bet=ud2(2)
      a=ud2(3)+ud2(4)*s
      ga=ud2(5)+ud2(6)*s**0.5
      gc=ud2(7)+ud2(8)*s
      b=ud2(9)+ud2(10)*s+ud2(11)*s**2
      gb=ud2(12)+ud2(13)*s+ud2(14)*s**2
      gd=ud2(15)+ud2(16)*s
      ge=ud2(17)+ud2(18)*s
      gep=ud2(19)+ud2(20)*s
      udpart2=suppr*grsf1(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
c
      alp=s2(1)
      bet=s2(2)
      a=s2(3)+s2(4)*s
      ga=s2(5)+s2(6)*s**0.5
      gc=s2(7)+s2(8)*s
      b=s2(9)+s2(10)*s+s2(11)*s**2
      gb=s2(12)+s2(13)*s+s2(14)*s**2
      gd=s2(15)+s2(16)*s
      ge=s2(17)+s2(18)*s
      gep=s2(19)+s2(20)*s
      spart2=suppr*grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
c
      alp=g2(1)
      bet=g2(2)
      a=g2(3)+g2(4)*s**0.5
      b=g2(5)+g2(6)*s**2
      gb=g2(7)+g2(8)*s
      ga=g2(9)+g2(10)*s**0.5+g2(11)*s**2
      gc=g2(12)+g2(13)*s**2
      gd=g2(14)+g2(15)*s+g2(16)*s**2
      ge=g2(17)+g2(18)*s
      gep=g2(19)+g2(20)*s
      gpart2=suppr*grsf1(x,s,alp,bet,a,b,ga,gb,gc,gd,ge,gep)
c
      ugam=upart1+udpart2
      DUPV = UGAM * ALPHEM
      DUSEA = DUPV
      dgam=dspart1+udpart2
      DDNV = DGAM * ALPHEM
      DDSEA = DDNV
      sgam=dspart1+spart2
      DSTR = SGAM * ALPHEM
      ggam=gpart1+gpart2
      DGL = GGAM * ALPHEM
C
      DCHM = 0.D0
      DBOT = 0.D0
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       FUNCTION GRVGF (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT REAL (A - Z)
       SX = SQRT (X)
       LX = ALOG (1./X)
       GRVGF  = (X**AK * (AG + BG * SX + C * X**BK)  +  S**AL
     1       * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       FUNCTION GRVGFS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT REAL (A - Z)
       IF (S .LE. SF) THEN
          GRVGFS = 0.0
       ELSE
          SX = SQRT (X)
          LX = ALOG (1./X)
          DS = S - SF
          GRVGFS = (DS * X**AK * (AG + BG * SX + C * X**BK) + DS**AL
     1         * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       END IF
       RETURN
       END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function grsf1(x,s,alp,bet,a,b,ga,gb,gc,gd,
     +                                ge,gep)
      implicit real*8 (a-h,o-z)
C
      grsf1=(x**a*(ga+gb*sqrt(x)+gc*x**b)+
     +      s**alp*exp(-ge+sqrt(gep*s**bet*log(1.d0/x))))*
     +      (1.d0-x)**gd
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,
     +                                ge,gep)
      implicit real*8 (a-h,o-z)
C
      grsf2=(s*x**a*(ga+gb*sqrt(x)+gc*x**b)+
     +      s**alp*exp(-ge+sqrt(gep*s**bet*log(1.d0/x))))*
     +      (1.d0-x)**gd
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


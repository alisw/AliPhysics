*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:27:41  mclareni
* Import version 8.04
*
* Revision 1.5  2000/04/17 18:54:47  plothow
* version 8.04: small mods, dlog to log
*
* Revision 1.4  2000/04/15 14:10:42  plothow
* Version 8.03: Nuclear PDFs introduced
*
* Revision 1.3  2000/04/09 10:06:18  plothow
* new version 8.00 - MRS 98 and 99, CTEQ5 and GRV98 added
*
* Revision 1.2  1996/10/28 16:01:42  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:50  plothow
* Version 7.01
*
*
*
* w50511c.inc
*
C
C------------------------------------------------------------------------------
C.
C     define the following sets by NPTYPE, NGROUP and NSET value :
C
C NPTYPE
C      NGROUP
C           NSET  Lambda   Name of set
C                 (Nfl=4)
C------------------------------------------------------------------------------
C             0    254  =  MRS Set (G) (MSbar)  (Default)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C   1    1    1         =  PRIVATE to user
C   1    1    2     45  =  BEBC                 [Z. Phys. C36 (1978) 1]
C   1    1    3    500  =  OR                   [Phys. Rev. D17 (1978) 3003]
C   1    1    4    470  =  BEP                  [Z. Phys. C2 (1979) 265]
C   1    1    5    400  =  GHR                  [Z. Phys. C13 (1982) 119]
C   1    1    6    200  =  DO Set 1             [Phys. Rev. D30 (1984) 49]
C   1    1    7    400  =  DO Set 2
C   1    1    8    200  =  EHLQ Set 1           [Rev.Mod.Phys. 56 (1984) 579,
C   1    1    9    290  =  EHLQ Set 2            Rev.Mod.Phys. 58 (1985) 1065]
C   1    1   10    177  =  New DO 1.1 Set LO    [Phys.Lett. 266B (1991) 126]
C------------------------------------------------------------------------------
C   1    2    1    200  =  DFLM soft valence    [Z. Phys. C39 (1988) 21]
C   1    2    2    200  =  DFLM hard valence
C   1    2    3    200  =  DFLM soft gluon
C   1    2    4    200  =  DFLM hard gluon
C   1    2    5    200  =  DFLM central average (LO)
C   1    2    6    300  =  DFLM central average (NLL - DIS)
C   1    2    7    160  =  DFLM 160    (DIS)    [unpublished]
C   1    2    8    260  =  DFLM 260    (DIS)
C   1    2    9    360  =  DFLM 360    (DIS)
C------------------------------------------------------------------------------
C   1    3    1    107  =  MRS Set 1            [Phys.Rev. D37 (1988) 1161]
C   1    3    2    250  =  MRS Set 2
C   1    3    3    178  =  MRS Set 3
C   1    3    4     91  =  MRS Set E            [Phys.Lett. 206B (1988) 327]
C   1    3    5    228  =  MRS Set B
C   1    3    6     91  =  MRS Set E'           [Mod.Phys.Lett. A4 (1989) 1135]
C   1    3    7    228  =  MRS Set B'
C   1    3    8    100  =  HMRS Set E  (Jan. 1990)   [retracted]
C   1    3    9    190  =  HMRS Set B  (Jan. 1990)
C   1    3   10    100  =  HMRS Set E  (March 1990)
C   1    3   11    190  =  HMRS Set B  (March 1990)
C   1    3   12    100  =  HMRS Set E+ (4.1990)  [Phys. Lett. 243B (1990) 421]
C   1    3   13    100  =  HMRS Set E- (4.1990)
C   1    3   14    100  =  HMRS Set E  (4.1990)  [Phys. Rev. D42 (1990) 798]
C   1    3   15    190  =  HMRS Set B  (4.1990)
C   1    3   16    100  =  HMRS Set B  (Aug. 1990)   [unpublished]
C   1    3   17    300  =  HMRS Set B  (Aug. 1990)
C   1    3   18    190  =  KMRS Set B- R2 Shad (7.1990) [Phys.Rev.D42(1990)3645
C   1    3   19    190  =  KMRS Set B- R5 Shad (7.1990)
C   1    3   20    190  =  KMRS Set B- (July 1990)
C   1    3   21    190  =  KMRS Set B0 (July 1990 - upgraded April 1992)
C   1    3   22    135  =  MRS Set B0-L135 (10.1990) [Phys. Rev D43 (1991) 3648
C   1    3   23    160  =  MRS Set B0-L160 (10.1990)
C   1    3   24    200  =  MRS Set B0-L200 (10.1990)
C   1    3   25    235  =  MRS Set B0-L235 (10.1990)
C   1    3   26    215  =  MRS Set S0 (4.1992) [retr.-Phys.Rev. D47 (1993) 867]
C   1    3   27    215  =  MRS Set D0 (4.1992) [retracted]
C   1    3   28    215  =  MRS Set D- (4.1992) [retracted]
C   1    3   29    230  =  MRS Set S0' (MSbar)      [Phys.Lett. 306B (1993) 145
C   1    3   30    230  =  MRS Set D0' (MSbar)
C   1    3   31    230  =  MRS Set D-' (MSbar)
C   1    3   32    230  =  MRS Set S0' (DIS)        [RAL preprint RAL-92-078]
C   1    3   33    230  =  MRS Set D0' (DIS)
C   1    3   34    230  =  MRS Set D-' (DIS)
C   1    3   35    230  =  MRS Set (H) (MSbar)
C   1    3   36    230  =  MRS Set (H) (DIS)
C   1    3   37    230  =  MRS Set (A) (MSbar)  [Phys.Rev. D50 (1994) 6734]
C   1    3   38    230  =  MRS Set (A) (MSbar) Fit
C   1    3   39    231  =  MRS Set (A')(MSbar)  [Phys.Lett. 354B (1995) 155]
C   1    3   40    231  =  MRS Set (A')(MSbar) Fit
C   1    3   41    254  =  MRS Set (G) (MSbar)  [RAL-95-021 (1995)]
C   1    3   42    254  =  MRS Set (G) (MSbar) Fit
C   1    3   43    230  =  MRS Set (A) (MSbar) low Q**2 [Phys.Rev. D51 (1995) 4756]
C   1    3   44    230  =  MRS Set (A) (DIS) low Q**2
C   1    3   45    150  =  MRS Set (A'-L150)(MSbar) [Phys.Lett. 356B (1995) 89]
C   1    3   46    201  =  MRS Set (A'-L201)(MSbar) 
C   1    3   47    266  =  MRS Set (A'-L266)(MSbar) 
C   1    3   48    344  =  MRS Set (A'-L344)(MSbar) 
C   1    3   49    435  =  MRS Set (A'-L435)(MSbar) 
C   1    3   50    542  =  MRS Set (A'-L542)(MSbar) 
C   1    3   51    344  =  MRS Set (J)  (MSbar)  [Phys.Lett. 381B (1996) 353]
C   1    3   52    507  =  MRS Set (J') (MSbar) 
C   1    3   53    241  =  MRS Set (R1) (MSbar)  [DTP/96/44 (1996)]
C   1    3   54    344  =  MRS Set (R2) (MSbar)  
C   1    3   55    241  =  MRS Set (R3) (MSbar)  
C   1    3   56    344  =  MRS Set (R4) (MSbar)  
C   1    3   57    241  =  MRS Set (R1) (MSbar) Fit
C   1    3   58    344  =  MRS Set (R2) (MSbar) Fit
C   1    3   59    300  =  MRRS Set (C1) (MSbar)  [DTP/96/102 (1996)]
C   1    3   60    300  =  MRRS Set (C2) (MSbar)
C   1    3   61    300  =  MRRS Set (C3) (MSbar)
C   1    3   62    300  =  MRST Set (c-g) (DIS) [DTP/98/10 (1998)]
C   1    3   63    300  =  MRST Set (h-g) (DIS) 
C   1    3   64    300  =  MRST Set (l-g) (DIS) 
C   1    3   65    229  =  MRST Set (l-a_s) (DIS) 
C   1    3   66    383  =  MRST Set (h-a_s) (DIS) 
C   1    3   67    300  =  MRST Set (c-g) (MSbar) [DTP/98/10 (1998)]
C   1    3   68    300  =  MRST Set (h-g) (MSbar) 
C   1    3   69    300  =  MRST Set (l-g) (MSbar) 
C   1    3   70    229  =  MRST Set (l-a_s) (MSbar) 
C   1    3   71    383  =  MRST Set (h-a_s) (MSbar) 
C   1    3   72    174  =  MRST Set (c-g) (lo) [DTP/98/52 (1998)]
C   1    3   73    174  =  MRST Set (h-g) (lo) 
C   1    3   74    174  =  MRST Set (l-g) (lo) 
C   1    3   75    136  =  MRST Set (l-a_s) (lo) 
C   1    3   76    216  =  MRST Set (h-a_s) (lo) 
C   1    3   77  292.5  =  MRST Set (ht) (MSbar) 
C   1    3   78    300  =  MRST Set (c-g) (DIS) [DTP/99/64 (1999)]
C   1    3   79    300  =  MRST Set (h-g) (DIS) 
C   1    3   80    300  =  MRST Set (l-g) (DIS) 
C   1    3   81    229  =  MRST Set (l-a_s) (DIS) 
C   1    3   82    383  =  MRST Set (h-a_s) (DIS) 
C   1    3   83  303.3  =  MRST Set (q-up)   (DIS) 
C   1    3   84  280.3  =  MRST Set (q-down) (DIS) 
C   1    3   85    300  =  MRST Set (s-up)   (DIS) 
C   1    3   86    300  =  MRST Set (s-down) (DIS) 
C   1    3   87    300  =  MRST Set (c-up)   (DIS) 
C   1    3   88    300  =  MRST Set (c-down) (DIS) 
C   1    3   89    300  =  MRST Set (c-g) (MSbar) 
C   1    3   90    300  =  MRST Set (h-g) (MSbar) 
C   1    3   91    300  =  MRST Set (l-g) (MSbar) 
C   1    3   92    229  =  MRST Set (l-a_s) (MSbar) 
C   1    3   93    383  =  MRST Set (h-a_s) (MSbar) 
C   1    3   94  303.3  =  MRST Set (q-up)   (MSbar) 
C   1    3   95  280.3  =  MRST Set (q-down) (MSbar) 
C   1    3   96    300  =  MRST Set (s-up)   (MSbar) 
C   1    3   97    300  =  MRST Set (s-down) (MSbar) 
C   1    3   98    300  =  MRST Set (c-up)   (MSbar) 
C   1    3   99    300  =  MRST Set (c-down) (MSbar) 
C   1    3  100    300  =  MRST Set (larger d/u) (MSbar) 
C------------------------------------------------------------------------------
C   1    4    1    212  =  MT set S1 (DIS)      [Z. Phys. C52 (1991) 13.
C   1    4    2    194  =  MT set B1 (DIS)
C   1    4    3    191  =  MT set B2 (DIS)
C   1    4    4    155  =  MT set E1 (DIS)
C   1    4    5    237  =  MT set 6  (1/2s) (DIS)
C   1    4    6    212  =  MT set S1 (MSbar)
C   1    4    7    194  =  MT set B1 (MSbar)
C   1    4    8    191  =  MT set B2 (MSbar)
C   1    4    9    155  =  MT set E1 (MSbar)
C   1    4   10    237  =  MT set 6  (1/2s) (MSbar)
C   1    4   11    144  =  MT set LO
C   1    4   12    168  =  CTEQ set 1LO      [retr.-Phys.Lett. 304B (1993) 159]
C   1    4   13    231  =  CTEQ set 1M  (MSbar)
C   1    4   14    231  =  CTEQ set 1MS (MSbar)
C   1    4   15    322  =  CTEQ set 1ML (MSbar)
C   1    4   16    247  =  CTEQ set 1D  (DIS)
C   1    4   17    190  =  CTEQ set 2LO           [retracted]
C   1    4   18    213  =  CTEQ set 2M  (MSbar)
C   1    4   19    208  =  CTEQ set 2MS (MSbar)
C   1    4   20    208  =  CTEQ set 2MF (MSbar)
C   1    4   21    322  =  CTEQ set 2ML (MSbar)
C   1    4   22    235  =  CTEQ set 2D  (DIS)
C   1    4   23    190  =  CTEQ set 2pL (LO)     [to be published]
C   1    4   24    213  =  CTEQ set 2pM  (MSbar)
C   1    4   25    208  =  CTEQ set 2pMS (MSbar)
C   1    4   26    208  =  CTEQ set 2pMF (MSbar)
C   1    4   27    322  =  CTEQ set 2pML (MSbar)
C   1    4   28    235  =  CTEQ set 2pD  (DIS)
C   1    4   29    177  =  CTEQ set 3LO          [MSUHEP-41024]
C   1    4   30    239  =  CTEQ set 3M   (MSbar)
C   1    4   31    247  =  CTEQ set 3D   (DIS)
C   1    4   32    236  =  CTEQ set 4L   (LO)    [Phys.Rev. D55 (1997) 1280]
C   1    4   33    298  =  CTEQ set 4D   (DIS)
C   1    4   34    298  =  CTEQ set 4M   (MSbar)
C   1    4   35    215  =  CTEQ set 4A1  (MSbar)
C   1    4   36    254  =  CTEQ set 4A2  (MSbar)
C   1    4   37    298  =  CTEQ set 4A3  (MSbar)
C   1    4   38    346  =  CTEQ set 4A4  (MSbar)
C   1    4   39    401  =  CTEQ set 4A5  (MSbar)
C   1    4   40    303  =  CTEQ set 4HJ  (MSbar)
C   1    4   41    261  =  CTEQ set 4lQ  (MSbar)
C   1    4   42    298  =  CTEQ set 4HQ  (MSbar)
C   1    4   43    298  =  CTEQ set 4HQ1 (MSbar)
C   1    4   44    385  =  CTEQ set 4F3  (NF=3)
C   1    4   45    292  =  CTEQ set 4F4  (NF=4)
C   1    4   46    192  =  CTEQ set 5L   (LO)    [hep-ph-9903282]
C   1    4   47    326  =  CTEQ set 5D   (DIS)
C   1    4   48    326  =  CTEQ set 5M   (MSbar)
C   1    4   49    326  =  CTEQ set 5HJ  (MSbar)
C   1    4   50    326  =  CTEQ set 5HQ  (MSbar)
C   1    4   51    395  =  CTEQ set 5F3  (NF=3)
C   1    4   52    309  =  CTEQ set 5F4  (NF=4)
C   1    4   53    326  =  CTEQ set 5M1  (MSbar)
C   1    4   54    326  =  CTEQ set 5HQ1 (MSbar)
C   1    4   55    192  =  CTEQ set 5Lp  (LO)   (parametrisation)
C   1    4   56    326  =  CTEQ set 5M1p (MSbar)(parametrisation)
C------------------------------------------------------------------------------
C   1    5    1    160  =  GRV Set HO           [retr.-Z. Phys. C48 (1990) 471]
C   1    5    2    220  =  GRV Set LO           [retracted]
C   1    5    3    200  =  GRV Set HO (MSbar)   [Z. Phys. C53 (1992) 127]
C   1    5    4    200  =  GRV Set LO
C   1    5    5    200  =  GRV Set LO           [Z. Phys. C67 (1995) 433]
C   1    5    6    200  =  GRV Set HO (MSbar)
C   1    5    7    200  =  GRV Set HO (DIS)
C   1    5    8    200  =  GRSV pol. LO stand   [Phys. Lett. 359B (1995) 201]
C   1    5    9    200  =  GRSV pol. LO val
C   1    5   10    200  =  GRSV pol. NLO stand (MSbar) [Phys. Rev. D53 (1996) 4775]
C   1    5   11    200  =  GRSV pol. NLO val   (MSbar)
C   1    5   12    200  =  GRV98 Set LO         [Eur. Phys. J. C5 (1998) 461]
C   1    5   13    200  =  GRV98 Set HO (MSbar)
C   1    5   14    200  =  GRV98 Set HO (DIS)
C------------------------------------------------------------------------------
C   1    6    1    230  =  ABFOW     (MSbar)    [Phys. Rev. D39 (1989) 3275]
C------------------------------------------------------------------------------
C   1    7    1    254  =  BM set A  (MSbar)    [Phys.Lett. 304B (1993) 318 and
C   1    7    2    254  =  BM set B  (MSbar)     CERN - TH 6739/92 (1992)]
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C   2    1    1    200  =  OW Set 1 (Pions)     [Phys. Rev. D30 (1984) 49]
C   2    1    2    400  =  OW Set 2 (Pions)
C------------------------------------------------------------------------------
C   2    3    1    190  =  SMRS-P Set 1 (Pions) [Phys. Rev. D45 (1992) 2349]
C   2    3    2    190  =  SMRS-P Set 2 (Pions)
C   2    3    3    190  =  SMRS-P Set 3 (Pions)
C------------------------------------------------------------------------------
C   2    5    1    200  =  GRV-P Set HO (Pions)  [Z. Phys. C53 (1992) 651]
C   2    5    2    200  =  GRV-P Set LO (Pions)
C------------------------------------------------------------------------------
C   2    6    1    231  =  ABFKW-P Set 1 (Pions) [Phys. Lett. 233B (1989) 517]
C   2    6    2    181  =  ABFKW-P Set 2 (Pions)
C   2    6    3    281  =  ABFKW-P Set 3 (Pions)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C   3    1    1    380  =  DO-G Set LO  (Photons) [Phys. Rev. D26 (1982) 1600]
C   3    1    2    440  =  DO-G Set HO  (Photons)
C------------------------------------------------------------------------------
C   3    2    1    400  =  DG-G Set 1 (Photons) [Z. Phys. C28 (1985) 451]
C   3    2    2    400  =  DG-G Set 2 (Photons)
C   3    2    3    400  =  DG-G Set 3 (Photons)
C   3    2    4    400  =  DG-G Set 4 (Photons)
C------------------------------------------------------------------------------
C   3    3    1    200  =  LAC-G Set 1 (Photons) [Phys.Lett. 269B (1991) 458]
C   3    3    2    200  =  LAC-G Set 2 (Photons)
C   3    3    3    200  =  LAC-G Set 3 (Photons)
C   3    3    3    200  =  GAL-G       (Photons) [hep-ph/9711355]
C------------------------------------------------------------------------------
C   3    4    1    200  =  GS-G Set HO  (Photons)   [Z. Phys. C56 (1992) 307]
C   3    4    2    200  =  GS-G LO Set 1 (Photons)
C   3    4    3    200  =  GS-G LO Set 2 (Photons)
C   3    4    4    200  =  GS-96-G HO (Photons)     [ANL-HEP-PR-96-33]
C   3    4    5    200  =  GS-96-G LO (Photons)
C------------------------------------------------------------------------------
C   3    5    1    200  =  GRV-G Set L-HO (Photons) [Phys. Rev. D46 (1992) 1973
C   3    5    2    200  =  GRV-G Set HO   (Photons) +Phys. Rev. D45 (1992) 3986
C   3    5    3    200  =  GRV-G Set LO   (Photons)
C   3    5    4    200  =  GRS-G Set LO   (Photons) [Phys. Rev. D51 (1995) 3220
C------------------------------------------------------------------------------
C   3    6    1    200  =  ACFGP-G Set HO (Photons)  [Z. Phys. C56 (1992) 589]
C   3    6    2    200  =  ACFGP-G Set HO-mc (Photons)
C   3    6    3    200  =  AFG-G Set HO (Photons)   [Z. Phys. C64 (1994) 621]
C------------------------------------------------------------------------------
C   3    8    1    200  =  WHIT-G Set 1 (Photons)  [Phys.Rev. D51 (1995) 3197]
C   3    8    2    200  =  WHIT-G Set 2 (Photons)
C   3    8    3    200  =  WHIT-G Set 3 (Photons)
C   3    8    4    200  =  WHIT-G Set 4 (Photons)
C   3    8    5    200  =  WHIT-G Set 5 (Photons)
C   3    8    6    200  =  WHIT-G Set 6 (Photons)
C------------------------------------------------------------------------------
C   3    9    1    200  =  SaS-G Set 1D (V1)       [Z. Phys. C68 (1995) 607]
C   3    9    2    200  =  SaS-G Set 1M (V1)
C   3    9    3    200  =  SaS-G Set 2D (V1)  
C   3    9    4    200  =  SaS-G Set 2M (V1)
C   3    9    5    200  =  SaS-G Set 1D (V2)       [Phys.Lett. 376B (1996) 193]
C   3    9    6    200  =  SaS-G Set 1M (V2)
C   3    9    7    200  =  SaS-G Set 2D (V2)  
C   3    9    8    200  =  SaS-G Set 2M (V2)
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C   4    1    1      0  =  EKS98 (Nuclear PDF's)   [Eur.Phys.J. C9 (1999) 61]
C------------------------------------------------------------------------------
C.

C-----------------------------------------------------------------------
#include "pdf/pilot.h"
C============================================================================
C                CTEQ Parton Distribution Functions: Version 4
C                          June 21, 1996
C
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C   Ref: MSUHEP-60426, CTEQ-604, e-Print Archive: hep-ph/9606399
C
C   This package contains 9 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C   Iset   PDF      Description             Alpha_s(Mz)  Q0(GeV)  Table_File
C ---------------------------------------------------------------------------
C   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
C   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
C   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
C   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
C   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
C   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
C   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
C   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
C   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
C   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ for which Q starts at a lower value of 0.7 GeV.  
C   The Table_Files are assumed to be in the working directory.
C   
C   The function Ctq4Fn (Iset, Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar)
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.

C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================
#if defined(CERNLIB_SINGLE)
      Function Ctq4Fn (Iset, Iparton, X, Q)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION FUNCTION Ctq4Fn (Iset, Iparton, X, Q)
#endif
#include "pdf/impdp.inc"
      Data Isetold, Isetmin, Isetmax / -987, 1, 10 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         if(Iset.eq.1)  Call Ctq4RM  (Iset)
         if(Iset.eq.2)  Call Ctq4RD  (Iset)
         if(Iset.eq.3)  Call Ctq4RL  (Iset)
         if(Iset.eq.4)  Call Ctq4RA1 (Iset)
         if(Iset.eq.5)  Call Ctq4RA2 (Iset)
         if(Iset.eq.6)  Call Ctq4RM  (Iset)
         if(Iset.eq.7)  Call Ctq4RA4 (Iset)
         if(Iset.eq.8)  Call Ctq4RA5 (Iset)
         if(Iset.eq.9)  Call Ctq4RHJ (Iset)
         if(Iset.eq.10) Call Ctq4RlQ (Iset)
         Isetold=Iset
      Endif
C
      Ctq4Fn = Ctq4PnX (Iparton, X, Q)
      if(Ctq4Fn.lt.0.D0)  Ctq4Fn = 0.D0
C
      Return
C                             ********************
      End

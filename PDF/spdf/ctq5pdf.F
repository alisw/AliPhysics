C------------------------------------------------------------------------------
#include "pdf/pilot.h"
C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.
C
C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------
C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================
#if defined(CERNLIB_SINGLE)
      Function Ctq5Pdf (Iset, Iparton, X, Q)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION Function Ctq5Pdf (Iset, Iparton, X, Q)
#endif
#include "pdf/impdp.inc"
      Data Isetold, Isetmin, Isetmax / -987, 1, 10 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         if(Iset.eq.1)  Call Ctq5M  (Iset)
         if(Iset.eq.2)  Call Ctq5D  (Iset)
         if(Iset.eq.3)  Call Ctq5L  (Iset)
         if(Iset.eq.4)  Call Ctq5HJ (Iset)
         if(Iset.eq.5)  Call Ctq5HQ (Iset)
         if(Iset.eq.6)  Call Ctq5F3 (Iset)
         if(Iset.eq.7)  Call Ctq5F4 (Iset)
         if(Iset.eq.8)  Call Ctq5M1 (Iset)
         if(Iset.eq.9)  Call Ctq5HQ1 (Iset)
         Isetold=Iset
      Endif
C
      Ctq5Pdf = Ctq5PnX (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0
C
      Return
C                             ********************
      End

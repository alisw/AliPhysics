*
* $Id$
*
* $Log$
* Revision 1.1.1.3  1997/04/15 14:02:00  cernlib
* Import version 7.08
*
* Revision 1.3  1997/04/08 15:53:25  pdflib
* Lambda_QCD for CTEQ4 sets corrected
*
* Revision 1.2  1996/10/28 15:03:30  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:25  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE PDFSET(PARM,VAL)
C
C **********************************************************************
C *                                                                    *
C *   Routine to set user parameters in the internal COMMON/W50511/,   *
C *                                                  COMMON/W50512/    *
C *                                              and COMMON/W50513/    *
C *   and provide filling of the internal            COMMON/W50518/    *
C *                                                  COMMON/W505110/   *
C *                                                  COMMON/W505120/   *
C *                                                  COMMON/W505121/   *
C *                                              and COMMON/W505122/   *
C *                                                                    *
C *   Author:  H. Plothow-Besch                                        *
C *            CERN-PPE, CH - 1211 Geneva 23, Switzerland              *
C *                                                                    *
C **********************************************************************
C
#include "pdf/impdp.inc"
C
#include "pdf/w5051p0.inc"
#include "pdf/w5051p1.inc"
#include "pdf/w5051p2.inc"
#include "pdf/w5051p3.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50513.inc"
#include "pdf/w50514.inc"
#include "pdf/w50514w.inc"
#include "pdf/w50515.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
#include "pdf/w50518.inc"
#include "pdf/w50519.inc"
#include "pdf/w505110.inc"
#include "pdf/w505120.inc"
#include "pdf/w505121.inc"
#include "pdf/w505122.inc"
#include "pdf/expdp.inc"
     +       PQCDL4(NPTYMX,NGRMAX,NSETMX), PQCDL5(NPTYMX,NGRMAX,NSETMX),
     +       PXMIN(NPTYMX,NGRMAX,NSETMX),  PXMAX(NPTYMX,NGRMAX,NSETMX),
     +       PQ2MIN(NPTYMX,NGRMAX,NSETMX), PQ2MAX(NPTYMX,NGRMAX,NSETMX)
C.
      PARAMETER (NSDIM = 16)
      CHARACTER*2 AV(NSDIM)
      CHARACTER*20 PARM(NCHDIM),STRING
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
C.
      EQUIVALENCE (VALUE,IVALUE)
C.
      DATA (PQCDL4(1,1,K),K=1,ND11) /
     +           0.0000D0, 0.0450D0, 0.5000D0, 0.4700D0, 0.4000D0,
     +           0.2000D0, 0.4000D0, 0.2000D0, 0.2900D0, 0.1770D0 /
      DATA (PQCDL4(1,2,K),K=1,ND12) /
     +           0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0,
     +           0.3000D0, 0.1600D0, 0.2600D0, 0.3600D0 /
      DATA (PQCDL4(1,3,K),K=1,ND13) /
     +           0.1070D0, 0.2500D0, 0.1780D0, 0.0910D0, 0.2280D0,
     +           0.0910D0, 0.2280D0, 0.1000D0, 0.1900D0, 0.1000D0,
     +           0.1900D0, 0.1000D0, 0.1000D0, 0.1000D0, 0.1900D0,
     +           0.1000D0, 0.3000D0, 0.1900D0, 0.1900D0, 0.1900D0,
     +           0.1900D0, 0.1350D0, 0.1600D0, 0.2000D0, 0.2350D0,
     +           0.2150D0, 0.2150D0, 0.2150D0, 0.2300D0, 0.2300D0,
     +           0.2300D0, 0.2300D0, 0.2300D0, 0.2300D0, 0.2300D0,
     +           0.2300D0, 0.2300D0, 0.2300D0, 0.2310D0, 0.2310D0,
     +           0.2540D0, 0.2540D0, 0.2300D0, 0.2300D0, 0.1500D0,
     +           0.2010D0, 0.2660D0, 0.3440D0, 0.4350D0, 0.5420D0,
     +           0.3440D0, 0.5070D0, 0.2410D0, 0.3440D0, 0.2410D0,
     +           0.3440D0, 0.2410D0, 0.3440D0 /
      DATA (PQCDL4(1,4,K),K=1,ND14) /
     +           0.2120D0, 0.1940D0, 0.1910D0, 0.1550D0, 0.2370D0,
     +           0.2120D0, 0.1940D0, 0.1910D0, 0.1550D0, 0.2370D0,
     +           0.1440D0, 0.1680D0, 0.2310D0, 0.2310D0, 0.3220D0,
     +           0.2470D0, 0.1900D0, 0.2130D0, 0.2080D0, 0.2080D0,
     +           0.3220D0, 0.2350D0, 0.1900D0, 0.2130D0, 0.2080D0,
     +           0.2080D0, 0.3220D0, 0.2350D0, 0.1770D0, 0.2390D0,
     +           0.2470D0, 0.2350D0, 0.2960D0, 0.2960D0, 0.2130D0,
     +           0.2530D0, 0.2960D0, 0.3440D0, 0.3990D0, 0.3020D0,
     +           0.2590D0 /
      DATA (PQCDL4(1,5,K),K=1,ND15) /
     +           0.1600D0, 0.2200D0, 0.2000D0, 0.2000D0, 0.2000D0,
     +           0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0,
     +           0.2000D0 /
      DATA (PQCDL4(1,6,K),K=1,ND16) /
     +           0.2300D0 /
      DATA (PQCDL4(1,7,K),K=1,ND17) /
     +           0.2540D0,0.2540D0 /
C.
      DATA (PQCDL4(2,1,K),K=1,ND21) /
     +           0.2000D0, 0.4000D0 /
      DATA (PQCDL4(2,3,K),K=1,ND23) /
     +           0.1900D0, 0.1900D0, 0.1900D0 /
      DATA (PQCDL4(2,5,K),K=1,ND25) /
     +           0.2000D0, 0.2000D0 /
      DATA (PQCDL4(2,6,K),K=1,ND26) /
     +           0.2310D0, 0.1810D0, 0.2810D0 /
C.
      DATA (PQCDL4(3,1,K),K=1,ND31) /
     +           0.3800D0, 0.4400D0 /
      DATA (PQCDL4(3,2,K),K=1,ND32) /
     +           0.4000D0, 0.4000D0, 0.4000D0, 0.4000D0 /
      DATA (PQCDL4(3,3,K),K=1,ND33) /
     +           0.2000D0, 0.2000D0, 0.2000D0 /
      DATA (PQCDL4(3,4,K),K=1,ND34) /
     +           0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0 /
      DATA (PQCDL4(3,5,K),K=1,ND35) /
     +           0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0 /
      DATA (PQCDL4(3,6,K),K=1,ND36) /
     +           0.2000D0, 0.2000D0, 0.2000D0 /
      DATA (PQCDL4(3,8,K),K=1,ND38) /
     +           0.4000D0, 0.4000D0, 0.4000D0, 0.4000D0, 0.4000D0,
     +           0.4000D0 /
      DATA (PQCDL4(3,9,K),K=1,ND39) /
     +           0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0, 0.2000D0,
     +           0.2000D0, 0.2000D0, 0.2000D0 /
C.
C.
      DATA (PQCDL5(1,1,K),K=1,ND11) /
     +           0.0000D0, 0.0300D0, 0.4111D0, 0.3844D0, 0.3226D0,
     +           0.1519D0, 0.3226D0, 0.1519D0, 0.2274D0, 0.1330D0 /
      DATA (PQCDL5(1,2,K),K=1,ND12) /
     +           0.1519D0, 0.1519D0, 0.1519D0, 0.1519D0, 0.1519D0,
     +           0.2048D0, 0.1016D0, 0.1745D0, 0.2513D0 /
      DATA (PQCDL5(1,3,K),K=1,ND13) /
     +           0.0651D0, 0.1671D0, 0.1144D0, 0.0544D0, 0.1507D0,
     +           0.0544D0, 0.1507D0, 0.0604D0, 0.1230D0, 0.0604D0,
     +           0.1230D0, 0.0604D0, 0.0604D0, 0.0604D0, 0.1230D0,
     +           0.0604D0, 0.2048D0, 0.1230D0, 0.1230D0, 0.1230D0,
     +           0.1230D0, 0.0842D0, 0.1016D0, 0.1303D0, 0.1559D0,
     +           0.1412D0, 0.1412D0, 0.1412D0, 0.1522D0, 0.1522D0,
     +           0.1522D0, 0.1522D0, 0.1522D0, 0.1522D0, 0.1522D0,
     +           0.1522D0, 0.1522D0, 0.1522D0, 0.1529D0, 0.1529D0,
     +           0.1700D0, 0.1700D0, 0.1522D0, 0.1522D0, 0.0946D0,
     +           0.1310D0, 0.1790D0, 0.2388D0, 0.3110D0, 0.3988D0,
     +           0.2388D0, 0.3698D0, 0.1604D0, 0.2388D0, 0.1604D0,
     +           0.2388D0, 0.1604D0, 0.2388D0 /
      DATA (PQCDL5(1,4,K),K=1,ND14) /
     +           0.1390D0, 0.1259D0, 0.1238D0, 0.0981D0, 0.1574D0,
     +           0.1390D0, 0.1259D0, 0.1238D0, 0.0981D0, 0.1574D0,
     +           0.1063D0, 0.1250D0, 0.1520D0, 0.1520D0, 0.2200D0,
     +           0.1640D0, 0.1430D0, 0.1390D0, 0.1350D0, 0.1350D0,
     +           0.2200D0, 0.1550D0, 0.1430D0, 0.1390D0, 0.1350D0,
     +           0.1350D0, 0.2200D0, 0.1550D0, 0.1320D0, 0.1580D0,
     +           0.1640D0, 0.1810D0, 0.2020D0, 0.2020D0, 0.1400D0,
     +           0.1690D0, 0.2020D0, 0.2390D0, 0.2820D0, 0.2060D0,
     +           0.1740D0 /
      DATA (PQCDL5(1,5,K),K=1,ND15) /
     +           0.1016D0, 0.1684D0, 0.1303D0, 0.1519D0, 0.1530D0,
     +           0.1310D0, 0.1310D0, 0.1530D0, 0.1530D0, 0.1310D0,
     +           0.1310D0 /
      DATA (PQCDL5(1,6,K),K=1,ND16) /
     +           0.1559D0 /
      DATA (PQCDL5(1,7,K),K=1,ND17) /
     +           0.1700D0,0.1700D0 /
C.
      DATA (PQCDL5(2,1,K),K=1,ND21) /
     +           0.1519D0, 0.3226D0 /
      DATA (PQCDL5(2,3,K),K=1,ND23) /
     +           0.1230D0, 0.1230D0, 0.1230D0 /
      DATA (PQCDL5(2,5,K),K=1,ND25) /
     +           0.1303D0, 0.1519D0 /
      DATA (PQCDL5(2,6,K),K=1,ND26) /
     +           0.1529D0, 0.1166D0, 0.1904D0 /
C.
      DATA (PQCDL5(3,1,K),K=1,ND31) /
     +           0.3051D0, 0.3150D0 /
      DATA (PQCDL5(3,2,K),K=1,ND32) /
     +           0.3226D0, 0.3226D0, 0.3226D0, 0.3226D0 /
      DATA (PQCDL5(3,3,K),K=1,ND33) /
     +           0.1519D0, 0.1519D0, 0.1519D0 /
      DATA (PQCDL5(3,4,K),K=1,ND34) /
     +           0.1303D0, 0.1519D0, 0.1519D0, 0.1303D0, 0.1519D0 /
      DATA (PQCDL5(3,5,K),K=1,ND35) /
     +           0.1303D0, 0.1303D0, 0.1519D0, 0.1530D0 /
      DATA (PQCDL5(3,6,K),K=1,ND36) /
     +           0.1303D0, 0.1303D0, 0.1303D0 /
      DATA (PQCDL5(3,8,K),K=1,ND38) /
     +           0.3226D0, 0.3226D0, 0.3226D0, 0.3226D0, 0.3226D0,
     +           0.3226D0 /
      DATA (PQCDL5(3,9,K),K=1,ND39) /
     +           0.3226D0, 0.3226D0, 0.3226D0, 0.3226D0, 0.3226D0,
     +           0.3226D0, 0.3226D0, 0.3226D0 /
C.
C.
      DATA (PXMIN(1,1,K),K=1,ND11) /
     +           0.0000D0, 5.000D-3, 5.000D-3, 5.000D-3, 1.000D-3,
     +           5.000D-3, 5.000D-3, 1.000D-4, 1.000D-4, 5.000D-5 /
      DATA (PXMIN(1,2,K),K=1,ND12) /
     +           5.000D-5, 5.000D-5, 5.000D-5, 5.000D-5, 5.000D-5,
     +           5.000D-5, 5.000D-5, 5.000D-5, 5.000D-5 /
      DATA (PXMIN(1,3,K),K=1,ND13) /
     +           1.000D-4, 1.000D-4, 1.000D-4, 1.000D-4, 1.000D-4,
     +           1.000D-4, 1.000D-4, 1.000D-4, 1.000D-4, 1.000D-4,
     +           1.000D-4, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5 /
      DATA (PXMIN(1,4,K),K=1,ND14) /
     +           1.00D-10, 1.00D-10, 1.00D-10, 1.00D-10, 1.00D-10,
     +           1.00D-10, 1.00D-10, 1.00D-10, 1.00D-10, 1.00D-10,
     +           1.00D-10, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5 /
      DATA (PXMIN(1,5,K),K=1,ND15) /
     +           1.000D-6, 1.000D-6, 1.000D-6, 1.000D-6, 1.000D-6,
     +           1.000D-6, 1.000D-6, 1.000D-4, 1.000D-4, 1.000D-4,
     +           1.000D-4 /
      DATA (PXMIN(1,6,K),K=1,ND16) /
     +           1.000D-4 /
      DATA (PXMIN(1,7,K),K=1,ND17) /
     +           1.000D-6, 1.000D-6 /
C.
      DATA (PXMIN(2,1,K),K=1,ND21) /
     +           5.000D-3, 5.000D-3 /
      DATA (PXMIN(2,3,K),K=1,ND23) /
     +           1.000D-5, 1.000D-5, 1.000D-5 /
      DATA (PXMIN(2,5,K),K=1,ND25) /
     +           1.000D-5, 1.000D-5 /
      DATA (PXMIN(2,6,K),K=1,ND26) /
     +           1.000D-3, 1.000D-3, 1.000D-3 /
C.
      DATA (PXMIN(3,1,K),K=1,ND31) /
     +           1.000D-5, 1.000D-5 /
      DATA (PXMIN(3,2,K),K=1,ND32) /
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5 /
      DATA (PXMIN(3,3,K),K=1,ND33) /
     +           1.000D-4, 1.000D-4, 1.000D-4 /
      DATA (PXMIN(3,4,K),K=1,ND34) /
     +           5.000D-4, 5.000D-4, 5.000D-4, 5.000D-4, 5.000D-4 /
      DATA (PXMIN(3,5,K),K=1,ND35) /
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5 /
      DATA (PXMIN(3,6,K),K=1,ND36) /
     +           1.370D-3, 1.370D-3, 1.370D-3 /
      DATA (PXMIN(3,8,K),K=1,ND38) /
     +           1.000D-3, 1.000D-3, 1.000D-3, 1.000D-3, 1.000D-3,
     +           1.000D-3 /
      DATA (PXMIN(3,9,K),K=1,ND39) /
     +           1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5, 1.000D-5,
     +           1.000D-5, 1.000D-5, 1.000D-5 /
C.
C.
      DATA (PXMAX(1,1,K),K=1,ND11) /
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99998D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0 /
      DATA (PXMAX(1,2,K),K=1,ND12) /
     +           0.95000D0, 0.95000D0, 0.95000D0, 0.95000D0, 0.95000D0,
     +           0.95000D0, 0.95000D0, 0.95000D0, 0.95000D0 /
      DATA (PXMAX(1,3,K),K=1,ND13) /
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0 /
      DATA (PXMAX(1,4,K),K=1,ND14) /
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0, 0.99998D0,
     +           0.99998D0 /
      DATA (PXMAX(1,5,K),K=1,ND15) /
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0, 0.99999D0,
     +           0.99999D0 /
      DATA (PXMAX(1,6,K),K=1,ND16) /
     +           0.99999D0 /
      DATA (PXMAX(1,7,K),K=1,ND17) /
     +           0.99999D0, 0.99999D0 /
C.
      DATA (PXMAX(2,1,K),K=1,ND21) /
     +           0.9998D0,  0.9998D0 /
      DATA (PXMAX(2,3,K),K=1,ND23) /
     +           0.9998D0,  0.9998D0,  0.9998D0 /
      DATA (PXMAX(2,5,K),K=1,ND25) /
     +           0.9998D0,  0.9998D0 /
      DATA (PXMAX(2,6,K),K=1,ND26) /
     +           0.9998D0,  0.9998D0,  0.9998D0 /
C.
      DATA (PXMAX(3,1,K),K=1,ND31) /
     +           0.9000D0, 0.9000D0 /
      DATA (PXMAX(3,2,K),K=1,ND32) /
     +           0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0 /
      DATA (PXMAX(3,3,K),K=1,ND33) /
     +           0.9998D0, 0.9998D0, 0.9998D0 /
      DATA (PXMAX(3,4,K),K=1,ND34) /
     +           0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0 /
      DATA (PXMAX(3,5,K),K=1,ND35) /
     +           0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0 /
      DATA (PXMAX(3,6,K),K=1,ND36) /
     +           0.9986D0, 0.9986D0, 0.9986D0 /
      DATA (PXMAX(3,8,K),K=1,ND38) /
     +           0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0,
     +           0.9998D0 /
      DATA (PXMAX(3,9,K),K=1,ND39) /
     +           0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0, 0.9998D0,
     +           0.9998D0, 0.9998D0, 0.9998D0 /
C.
C.
      DATA (PQ2MIN(1,1,K),K=1,ND11) /
     +           0.000D0, 0.500D0, 1.800D0, 5.000D0, 4.000D0,
     +           4.000D0, 4.000D0, 5.000D0, 5.000D0, 4.000D0 /
      DATA (PQ2MIN(1,2,K),K=1,ND12) /
     +           10.00D0, 10.00D0, 10.00D0, 10.00D0, 10.00D0,
     +           10.00D0, 10.00D0, 10.00D0, 10.00D0 /
      DATA (PQ2MIN(1,3,K),K=1,ND13) /
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 0.625D0, 0.625D0, 5.000D0,
     +           5.000D0, 5.000D0, 5.000D0, 5.000D0, 5.000D0,
     +           5.000D0, 5.000D0, 1.250D0, 1.250D0, 1.250D0,
     +           1.250D0, 2.000D0, 2.000D0 /
      DATA (PQ2MIN(1,4,K),K=1,ND14) /
     +           4.000D0, 4.000D0, 4.000D0, 4.000D0, 4.000D0,
     +           4.000D0, 4.000D0, 4.000D0, 4.000D0, 4.000D0,
     +           4.000D0, 4.000D0, 4.000D0, 4.000D0, 4.000D0,
     +           4.000D0, 4.000D0, 4.000D0, 4.000D0, 4.000D0,
     +           4.000D0, 4.000D0, 4.000D0, 4.000D0, 4.000D0,
     +           4.000D0, 4.000D0, 4.000D0, 4.000D0, 4.000D0,
     +           4.000D0, 2.560D0, 2.560D0, 2.560D0, 2.560D0,
     +           2.560D0, 2.560D0, 2.560D0, 2.560D0, 2.560D0,
     +           0.490D0 /
      DATA (PQ2MIN(1,5,K),K=1,ND15) /
     +           0.200D0, 0.200D0, 0.300D0, 0.250D0, 0.400D0,
     +           0.400D0, 0.400D0, 0.400D0, 0.400D0, 0.400D0,
     +           0.400D0 /
      DATA (PQ2MIN(1,6,K),K=1,ND16) /
     +           2.0000D0 /
      DATA (PQ2MIN(1,7,K),K=1,ND17) /
     +           5.000D0, 5.000D0 /
C.
      DATA (PQ2MIN(2,1,K),K=1,ND21) /
     +           4.0000D0, 4.0000D0 /
      DATA (PQ2MIN(2,3,K),K=1,ND23) /
     +           5.0000D0, 5.0000D0, 5.0000D0 /
      DATA (PQ2MIN(2,5,K),K=1,ND25) /
     +           0.3000D0, 0.2500D0 /
      DATA (PQ2MIN(2,6,K),K=1,ND26) /
     +           2.0000D0, 2.0000D0, 2.0000D0 /
C.
      DATA (PQ2MIN(3,1,K),K=1,ND31) /
     +           10.000D0, 10.000D0 /
      DATA (PQ2MIN(3,2,K),K=1,ND32) /
     +           1.0000D0, 1.0000D0, 20.000D0, 200.00D0 /
      DATA (PQ2MIN(3,3,K),K=1,ND33) /
     +           4.0000D0, 4.0000D0, 1.0000D0 /
      DATA (PQ2MIN(3,4,K),K=1,ND34) /
     +           5.3000D0, 5.3000D0, 5.3000D0, 5.3000D0, 5.3000D0 /
      DATA (PQ2MIN(3,5,K),K=1,ND35) /
     +           0.3000D0, 0.3000D0, 0.2500D0, 0.6000D0 /
      DATA (PQ2MIN(3,6,K),K=1,ND36) /
     +           2.0000D0, 2.0000D0, 2.0000D0 /
      DATA (PQ2MIN(3,8,K),K=1,ND38) /
     +           4.0000D0, 4.0000D0, 4.0000D0, 4.0000D0, 4.0000D0,
     +           4.0000D0 /
      DATA (PQ2MIN(3,9,K),K=1,ND39) /
     +           0.3600D0, 0.3600D0, 4.0000D0, 4.0000D0, 0.3600D0,
     +           0.3600D0, 4.0000D0, 4.0000D0 /
C.
C.
      DATA (PQ2MAX(1,1,K),K=1,ND11) /
     +           0.000D0, 1.000D6, 1.000D6, 1.000D6, 1.000D6,
     +           1.000D6, 1.000D6, 1.000D8, 1.000D8, 4.000D7 /
      DATA (PQ2MAX(1,2,K),K=1,ND12) /
     +           1.000D8, 1.000D8, 1.000D8, 1.000D8, 1.000D8,
     +           1.000D8, 1.000D8, 1.000D8, 1.000D8 /
      DATA (PQ2MAX(1,3,K),K=1,ND13) /
     +           1.310D6, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           1.310D6, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           1.310D6, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           1.310D6, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           2.680D9, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           2.680D9, 2.680D9, 2.680D9, 2.680D9, 2.680D9,
     +           2.680D9, 2.680D9, 2.680D9, 2.680D9, 2.680D9,
     +           2.680D9, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           1.310D6, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           1.310D6, 1.310D6, 1.310D6, 1.310D6, 1.310D6,
     +           1.310D6, 1.310D6, 1.000D7, 1.000D7, 1.000D7,
     +           1.000D7, 1.000D7, 1.000D7 /
      DATA (PQ2MAX(1,4,K),K=1,ND14) /
     +           1.000D8, 1.000D8, 1.000D8, 1.000D8, 1.000D8,
     +           1.000D8, 1.000D8, 1.000D8, 1.000D8, 1.000D8,
     +           1.000D8, 1.000D6, 1.000D6, 1.000D6, 1.000D6,
     +           1.000D6, 1.000D6, 1.000D6, 1.000D6, 1.000D6,
     +           1.000D6, 1.000D6, 1.000D6, 1.000D6, 1.000D6,
     +           1.000D6, 1.000D6, 1.000D6, 1.000D6, 1.000D6,
     +           1.000D6, 1.000D8, 1.000D8, 1.000D8, 1.000D8,
     +           1.000D6, 1.000D8, 1.000D8, 1.000D8, 1.000D8,
     +           1.000D8 /
      DATA (PQ2MAX(1,5,K),K=1,ND15) /
     +           1.000D8, 1.000D8, 1.000D8, 1.000D8, 1.000D6,
     +           1.000D6, 1.000D6, 1.000D4, 1.000D4, 1.000D4,
     +           1.000D4 /
      DATA (PQ2MAX(1,6,K),K=1,ND16) /
     +           1.0000D8 /
      DATA (PQ2MAX(1,7,K),K=1,ND17) /
     +           1.310D6, 1.310D6 /
C.
      DATA (PQ2MAX(2,1,K),K=1,ND21) /
     +           2.0000D3, 2.0000D3 /
      DATA (PQ2MAX(2,3,K),K=1,ND23) /
     +           1.3100D6, 1.3100D6, 1.3100D6 /
      DATA (PQ2MAX(2,5,K),K=1,ND25) /
     +           1.0000D6, 1.0000D6 /
      DATA (PQ2MAX(2,6,K),K=1,ND26) /
     +           1.0000D8, 1.0000D8, 1.0000D8 /
C.
      DATA (PQ2MAX(3,1,K),K=1,ND31) /
     +           1.0000D4, 1.0000D4 /
      DATA (PQ2MAX(3,2,K),K=1,ND32) /
     +           1.0000D4, 50.000D0, 500.00D0, 1.0000D4 /
      DATA (PQ2MAX(3,3,K),K=1,ND33) /
     +           1.0000D5, 1.0000D5, 1.0000D5 /
      DATA (PQ2MAX(3,4,K),K=1,ND34) /
     +           1.0000D8, 1.0000D8, 1.0000D8, 1.0000D8, 1.0000D8 /
      DATA (PQ2MAX(3,5,K),K=1,ND35) /
     +           1.0000D6, 1.0000D6, 1.0000D6, 5.0000D4 /
      DATA (PQ2MAX(3,6,K),K=1,ND36) /
     +           5.5000D5, 5.5000D5, 5.5000D5/
      DATA (PQ2MAX(3,8,K),K=1,ND38) /
     +           2.5000D3, 2.5000D3, 2.5000D3, 2.5000D3, 2.5000D3,
     +           2.5000D3 /
      DATA (PQ2MAX(3,9,K),K=1,ND39) /
     +           5.0000D4, 5.0000D4, 5.0000D4, 5.0000D4, 5.0000D4,
     +           5.0000D4, 5.0000D4, 5.0000D4 /
C.
      DATA AV /'BE','OR','GH','DO','EH','DF','MR','HM','KM',
     +         'MT','CT','GR','AB','BM','WI','SA'/
C.
      DATA ZEROD/0.D0/,PONED/0.1D0/,ONED/1.D0/,ONEDO/1.1D0/,TWOD/2.D0/
      DATA ZEROMIN/0.0001D0/
      DATA ISTART/0/
C.
      SAVE ISTART, PQCDL4, PQCDL5, PXMIN, PXMAX, PQ2MIN, PQ2MAX,
     +     /W50511/, /W505110/, /W505120/, /W505121/, /W505122/,
     +     /W50514/, /W50514W/, /W50515/, /W50516/, /W50518/
C.
      IF (ISTART.EQ.0) IFLSET = 0
C.
C Set the PDFLIB parameters to the values of the previous call to PDFSET
      IF(IFLSET.EQ.1) THEN
         IFLPRT = IFLPRTP
         NPTYPE = NPTYPP
         NGROUP = NGROPP
         NSET   = NSETP
         NFL    = NFLP
         LO     = LOP
         TMAS   = TMASP
         QCDL4  = QCDL4P
         QCDL5  = QCDL5P
         XMIN   = XMINP
         XMAX   = XMAXP
         Q2MIN  = Q2MINP
         Q2MAX  = Q2MAXP
      ENDIF
C.
C Initialisation of PDFLIB parameters (first call to PDFSET)
      IF(ISTART.EQ.0) THEN
C.
C set maximum dimensions for PDFs per particle type and per author group
         CALL UFILL(NPGSMX,1,NPTYMX*NGRMAX,0)
         NPGSMX(1,1) = ND11
         NPGSMX(1,2) = ND12
         NPGSMX(1,3) = ND13
         NPGSMX(1,4) = ND14
         NPGSMX(1,5) = ND15
         NPGSMX(1,6) = ND16
         NPGSMX(1,7) = ND17
         NPGSMX(1,8) = ND18
         NPGSMX(1,9) = ND19
         NPGSMX(2,1) = ND21
         NPGSMX(2,2) = ND22
         NPGSMX(2,3) = ND23
         NPGSMX(2,4) = ND24
         NPGSMX(2,5) = ND25
         NPGSMX(2,6) = ND26
         NPGSMX(2,7) = ND27
         NPGSMX(2,8) = ND28
         NPGSMX(2,9) = ND29
         NPGSMX(3,1) = ND31
         NPGSMX(3,2) = ND32
         NPGSMX(3,3) = ND33
         NPGSMX(3,4) = ND34
         NPGSMX(3,5) = ND35
         NPGSMX(3,6) = ND36
         NPGSMX(3,7) = ND37
         NPGSMX(3,8) = ND38
         NPGSMX(3,9) = ND39
C.
C set retracted PDF sets invalid
         CALL UFILL(NSETFL,1,NPTYMX*NGRMAX*NSETMX,0)
         DO KPTY = 1,NPTYMX
            DO KGRO = 1,NGRMAX
               IF(NPGSMX(KPTY,KGRO) .EQ. 0) GOTO 2
               DO ISET = 1,NPGSMX(KPTY,KGRO)
                  NSETFL(KPTY,KGRO,ISET) = 1
                  IF (KPTY.EQ.1) THEN
                     IF (KGRO.EQ.3) THEN
                        IF (ISET.EQ.8) THEN
                           NSETFL(KPTY,KGRO,ISET) = 0
                        ENDIF
                     ELSEIF (KGRO.EQ.5) THEN
                        IF (ISET.EQ.1 .OR. ISET.EQ.2) THEN
                           NSETFL(KPTY,KGRO,ISET) = 0
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
    2          CONTINUE
            ENDDO
         ENDDO
C.
C set character*8 NAME for each PDF
         DO NPTY=1,NPTYMX
            DO NGR=1,NGRMAX
               DO NSET=1,NSETMX
                  SFNAME(NPTY,NGR,NSET) = '        '
               ENDDO
            ENDDO
         ENDDO
         SFNAME(1,1,1)  = 'PRIVATE '
         SFNAME(1,1,2)  = 'BEBC    '
         SFNAME(1,1,3)  = 'OR      '
         SFNAME(1,1,4)  = 'BEP     '
         SFNAME(1,1,5)  = 'GHR     '
         SFNAME(1,1,6)  = 'DO-1    '
         SFNAME(1,1,7)  = 'DO-2    '
         SFNAME(1,1,8)  = 'EHLQ-1  '
         SFNAME(1,1,9)  = 'EHLQ-2  '
         SFNAME(1,1,10) = 'DO-1.1  '
         SFNAME(1,2,1)  = 'DFLM-SV '
         SFNAME(1,2,2)  = 'DFLM-HV '
         SFNAME(1,2,3)  = 'DFLM-SG '
         SFNAME(1,2,4)  = 'DFLM-HG '
         SFNAME(1,2,5)  = 'DFLM-LO '
         SFNAME(1,2,6)  = 'DFLM-HO '
         SFNAME(1,2,7)  = 'DFLM-160'
         SFNAME(1,2,8)  = 'DFLM-260'
         SFNAME(1,2,9)  = 'DFLM-360'
         SFNAME(1,3,1)  = 'MRS-1   '
         SFNAME(1,3,2)  = 'MRS-2   '
         SFNAME(1,3,3)  = 'MRS-3   '
         SFNAME(1,3,4)  = 'MRS-E   '
         SFNAME(1,3,5)  = 'MRS-B   '
         SFNAME(1,3,6)  = 'MRS-EP  '
         SFNAME(1,3,7)  = 'MRS-BP  '
         SFNAME(1,3,8)  = 'HMRS-E1.'
         SFNAME(1,3,9)  = 'HMRS-B1.'
         SFNAME(1,3,10) = 'HMRS-E3.'
         SFNAME(1,3,11) = 'HMRS-B3.'
         SFNAME(1,3,12) = 'HMRSE+  '
         SFNAME(1,3,13) = 'HMRSE-  '
         SFNAME(1,3,14) = 'HMRS-E  '
         SFNAME(1,3,15) = 'HMRS-B  '
         SFNAME(1,3,16) = 'HMRSB100'
         SFNAME(1,3,17) = 'HMRSB300'
         SFNAME(1,3,18) = 'KMRS-BM2'
         SFNAME(1,3,19) = 'KMRS-BM5'
         SFNAME(1,3,20) = 'KMRS-BM '
         SFNAME(1,3,21) = 'KMRS-B0 '
         SFNAME(1,3,22) = 'MRS-B0-1'
         SFNAME(1,3,23) = 'MRS-B0-2'
         SFNAME(1,3,24) = 'MRS-B0-3'
         SFNAME(1,3,25) = 'MRS-B0-4'
         SFNAME(1,3,26) = 'MRS-S0  '
         SFNAME(1,3,27) = 'MRS-D0  '
         SFNAME(1,3,28) = 'MRS-DM  '
         SFNAME(1,3,29) = 'MRS-S0P '
         SFNAME(1,3,30) = 'MRS-D0P '
         SFNAME(1,3,31) = 'MRS-DMP '
         SFNAME(1,3,32) = 'MRS-S0PD'
         SFNAME(1,3,33) = 'MRS-D0PD'
         SFNAME(1,3,34) = 'MRS-DMPD'
         SFNAME(1,3,35) = 'MRS-(H) '
         SFNAME(1,3,36) = 'MRS-(H)D'
         SFNAME(1,3,37) = 'MRS-(A) '
         SFNAME(1,3,38) = 'MRS-(A)F'
         SFNAME(1,3,39) = 'MRS-Ap  '
         SFNAME(1,3,40) = 'MRS-Ap F'
         SFNAME(1,3,41) = 'MRS-G   '
         SFNAME(1,3,42) = 'MRS-GFit'
         SFNAME(1,3,43) = 'MRS-A lQ'
         SFNAME(1,3,44) = 'MRS-ADlQ'
         SFNAME(1,3,45) = 'MRS-Ap-1'
         SFNAME(1,3,46) = 'MRS-Ap-2'
         SFNAME(1,3,47) = 'MRS-Ap-3'
         SFNAME(1,3,48) = 'MRS-Ap-4'
         SFNAME(1,3,49) = 'MRS-Ap-5'
         SFNAME(1,3,50) = 'MRS-Ap-6'
         SFNAME(1,3,51) = 'MRS-J   '
         SFNAME(1,3,52) = 'MRS-Jp  '
         SFNAME(1,3,53) = 'MRS-R1  '
         SFNAME(1,3,54) = 'MRS-R2  '
         SFNAME(1,3,55) = 'MRS-R3  '
         SFNAME(1,3,56) = 'MRS-R4  '
         SFNAME(1,3,57) = 'MRS-R1Ft'
         SFNAME(1,3,58) = 'MRS-R2Ft'
         SFNAME(1,4,1)  = 'MT-S1DIS'
         SFNAME(1,4,2)  = 'MT-B1DIS'
         SFNAME(1,4,3)  = 'MT-B2DIS'
         SFNAME(1,4,4)  = 'MT-E1DIS'
         SFNAME(1,4,5)  = 'MT-6 DIS'
         SFNAME(1,4,6)  = 'MT-S1MSB'
         SFNAME(1,4,7)  = 'MT-B1MSB'
         SFNAME(1,4,8)  = 'MT-B2MSB'
         SFNAME(1,4,9)  = 'MT-E1MSB'
         SFNAME(1,4,10) = 'MT-6 MSB'
         SFNAME(1,4,11) = 'MT-LO   '
         SFNAME(1,4,12) = 'CTEQ-1LO'
         SFNAME(1,4,13) = 'CTEQ-1M '
         SFNAME(1,4,14) = 'CTEQ-1MS'
         SFNAME(1,4,15) = 'CTEQ-1ML'
         SFNAME(1,4,16) = 'CTEQ-1D '
         SFNAME(1,4,17) = 'CTEQ-2LO'
         SFNAME(1,4,18) = 'CTEQ-2M '
         SFNAME(1,4,19) = 'CTEQ-2MS'
         SFNAME(1,4,20) = 'CTEQ-2MF'
         SFNAME(1,4,21) = 'CTEQ-2ML'
         SFNAME(1,4,22) = 'CTEQ-2D '
         SFNAME(1,4,23) = 'CTEQ-2pL'
         SFNAME(1,4,24) = 'CTEQ-2pM'
         SFNAME(1,4,25) = 'CTEQ2pMS'
         SFNAME(1,4,26) = 'CTEQ2pMF'
         SFNAME(1,4,27) = 'CTEQ2pML'
         SFNAME(1,4,28) = 'CTEQ-2pD'
         SFNAME(1,4,29) = 'CTEQ-3L '
         SFNAME(1,4,30) = 'CTEQ-3M '
         SFNAME(1,4,31) = 'CTEQ-3D '
         SFNAME(1,4,32) = 'CTEQ-4L '
         SFNAME(1,4,33) = 'CTEQ-4D '
         SFNAME(1,4,34) = 'CTEQ-4M '
         SFNAME(1,4,35) = 'CTEQ-4A1'
         SFNAME(1,4,36) = 'CTEQ-4A2'
         SFNAME(1,4,37) = 'CTEQ-4A3'
         SFNAME(1,4,38) = 'CTEQ-4A4'
         SFNAME(1,4,39) = 'CTEQ-4A5'
         SFNAME(1,4,40) = 'CTEQ-4HJ'
         SFNAME(1,4,41) = 'CTEQ-4lQ'
         SFNAME(1,5,1)  = 'GRV-HO  '
         SFNAME(1,5,2)  = 'GRV-LO  '
         SFNAME(1,5,3)  = 'GRV-HO  '
         SFNAME(1,5,4)  = 'GRV-LO  '
         SFNAME(1,5,5)  = 'GRV94-LO'
         SFNAME(1,5,6)  = 'GRV94-MS'
         SFNAME(1,5,7)  = 'GRV94-DI'
         SFNAME(1,5,8)  = 'GRSV-LOS'
         SFNAME(1,5,9)  = 'GRSV-LOV'
         SFNAME(1,5,10) = 'GRSV-HOS'
         SFNAME(1,5,11) = 'GRSV-HOV'
         SFNAME(1,6,1)  = 'ABFOW   '
         SFNAME(1,7,1)  = 'BM-A    '
         SFNAME(1,7,2)  = 'BM-B    '
         SFNAME(2,1,1)  = 'OW-P-1  '
         SFNAME(2,1,2)  = 'OW-P-2  '
         SFNAME(2,3,1)  = 'SMRS-P-1'
         SFNAME(2,3,2)  = 'SMRS-P-2'
         SFNAME(2,3,3)  = 'SMRS-P-3'
         SFNAME(2,5,1)  = 'GRV-P-HO'
         SFNAME(2,5,2)  = 'GRV-P-LO'
         SFNAME(2,6,1)  = 'ABFKW-P1'
         SFNAME(2,6,2)  = 'ABFKW-P2'
         SFNAME(2,6,3)  = 'ABFKW-P3'
         SFNAME(3,1,1)  = 'DO-G-LO '
         SFNAME(3,1,2)  = 'DO-G-HO '
         SFNAME(3,2,1)  = 'DG-G-1  '
         SFNAME(3,2,2)  = 'DG-G-2  '
         SFNAME(3,2,3)  = 'DG-G-3  '
         SFNAME(3,2,4)  = 'DG-G-4  '
         SFNAME(3,3,1)  = 'LAC-G-1 '
         SFNAME(3,3,2)  = 'LAC-G-2 '
         SFNAME(3,3,3)  = 'LAC-G-3 '
         SFNAME(3,4,1)  = 'GS-G-1  '
         SFNAME(3,4,2)  = 'GS-G-2  '
         SFNAME(3,4,3)  = 'GS-G-3  '
         SFNAME(3,4,4)  = 'GS-96-HO'
         SFNAME(3,4,5)  = 'GS-96-LO'
         SFNAME(3,5,1)  = 'GRV-GLHO'
         SFNAME(3,5,2)  = 'GRV-G-HO'
         SFNAME(3,5,3)  = 'GRV-G-LO'
         SFNAME(3,5,4)  = 'GRS-G-LO'
         SFNAME(3,6,1)  = 'ACFGP-GH'
         SFNAME(3,6,2)  = 'ACFGP-GC'
         SFNAME(3,6,3)  = 'AFG-G-HO'
         SFNAME(3,8,1)  = 'WHIT-G-1'
         SFNAME(3,8,2)  = 'WHIT-G-2'
         SFNAME(3,8,3)  = 'WHIT-G-3'
         SFNAME(3,8,4)  = 'WHIT-G-4'
         SFNAME(3,8,5)  = 'WHIT-G-5'
         SFNAME(3,8,6)  = 'WHIT-G-6'
         SFNAME(3,9,1)  = 'SaS-1G1D'
         SFNAME(3,9,2)  = 'SaS-1G1M'
         SFNAME(3,9,3)  = 'SaS-1G2D'
         SFNAME(3,9,4)  = 'SaS-1G2M'
         SFNAME(3,9,5)  = 'SaS-2G1D'
         SFNAME(3,9,6)  = 'SaS-2G1M'
         SFNAME(3,9,7)  = 'SaS-2G2D'
         SFNAME(3,9,8)  = 'SaS-2G2M'
C.
C set cross reference tables for NEW PDFLIB format version 4.0
C       (give OLD to get NEW)
         DO I=1,MODEMX
         NPTYCR(I) = I/100 + 1
         IF(I.LE.99) THEN
           IF(I.LE.  9) THEN
              NGROCR(I) = 1
              NSETCR(I) = I
           ELSEIF(I.GE. 10 .AND. I.LE. 19) THEN
              NGROCR(I) = 2
              NSETCR(I) = I - 9
           ELSEIF(I.GE. 20 .AND. I.LE. 47) THEN
              NGROCR(I) = 3
              NSETCR(I) = I - 19
           ELSEIF(I.GE. 48 .AND. I.LE. 49) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE. 50 .AND. I.LE. 61) THEN
              NGROCR(I) = 4
              NSETCR(I) = I - 49
           ELSEIF(I.GE. 62 .AND. I.LE. 69) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE. 70 .AND. I.LE. 73) THEN
              NGROCR(I) = 5
              NSETCR(I) = I - 69
           ELSEIF(I.EQ. 74 )               THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.EQ. 75 )               THEN
              NGROCR(I) = 1
              NSETCR(I) = I - 65
           ELSEIF(I.GE. 76 .AND. I.LE. 79) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.EQ. 80 )               THEN
              NGROCR(I) = 6
              NSETCR(I) = I - 79
           ELSEIF(I.GE. 81 .AND. I.LE. 85) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE. 86 .AND. I.LE. 87) THEN
              NGROCR(I) = 7
              NSETCR(I) = I - 85
           ELSEIF(I.GE. 88 .AND. I.LE. 99) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ENDIF
         ELSEIF(I.GE.100 .AND. I.LE.199) THEN
           IF(I.GE.100 .AND. I.LE.105) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.106 .AND. I.LE.107) THEN
              NGROCR(I) = 1
              NSETCR(I) = I - 105
           ELSEIF(I.GE.108 .AND. I.LE.129) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.130 .AND. I.LE.132) THEN
              NGROCR(I) = 3
              NSETCR(I) = I - 129
           ELSEIF(I.GE.133 .AND. I.LE.171) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.172 .AND. I.LE.173) THEN
              NGROCR(I) = 5
              NSETCR(I) = I - 171
           ELSEIF(I.GE.174 .AND. I.LE.179) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.180 .AND. I.LE.182) THEN
              NGROCR(I) = 6
              NSETCR(I) = I - 179
           ELSEIF(I.GE.183 .AND. I.LE.199) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ENDIF
         ELSEIF(I.GE.200 .AND. I.LE.(MODEMX-1)) THEN
           IF(I.GE.200 .AND. I.LE.205) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.206 .AND. I.LE.207) THEN
              NGROCR(I) = 1
              NSETCR(I) = I - 205
           ELSEIF(I.GE.208 .AND. I.LE.209) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.210 .AND. I.LE.213) THEN
              NGROCR(I) = 2
              NSETCR(I) = I - 209
           ELSEIF(I.GE.214 .AND. I.LE.219) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.220 .AND. I.LE.222) THEN
              NGROCR(I) = 3
              NSETCR(I) = I - 219
           ELSEIF(I.GE.223 .AND. I.LE.229) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.230 .AND. I.LE.232) THEN
              NGROCR(I) = 4
              NSETCR(I) = I - 229
           ELSEIF(I.GE.233 .AND. I.LE.270) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.271 .AND. I.LE.273) THEN
              NGROCR(I) = 5
              NSETCR(I) = I - 270
           ELSEIF(I.GE.274 .AND. I.LE.279) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ELSEIF(I.GE.280 .AND. I.LE.281) THEN
              NGROCR(I) = 6
              NSETCR(I) = I - 279
           ELSEIF(I.GE.282 .AND. I.LE.(MODEMX-1)) THEN
              NGROCR(I) = -1
              NSETCR(I) = -1
           ENDIF
         ELSEIF(I.EQ.MODEMX) THEN
           NPTYCR(I) = -1
           NGROCR(I) = -1
           NSETCR(I) = -1
         ENDIF
      ENDDO
C.
C Now get opposite cross reference
C       (give NEW to get OLD)
         CALL UFILL(MODECR,1,NPTYMX*NGRMAX*NSETMX,-1)
         DO I = 1,NPTYMX
            DO J = 1,NGRMAX
               IF (I.EQ.1) THEN
                  IF (J.EQ.1) KT =  0
                  IF (J.EQ.2) KT =  9
                  IF (J.EQ.3) KT = 19
                  IF (J.EQ.4) KT = 49
                  IF (J.EQ.5) KT = 69
                  IF (J.EQ.6) KT = 79
                  IF (J.EQ.7) KT = 85
               ELSEIF (I.EQ.2) THEN
                  IF (J.EQ.1) KT = 105
                  IF (J.EQ.3) KT = 129
                  IF (J.EQ.5) KT = 171
                  IF (J.EQ.6) KT = 179
               ELSEIF (I.EQ.3) THEN
                  IF (J.EQ.1) KT = 205
                  IF (J.EQ.2) KT = 209
                  IF (J.EQ.3) KT = 219
                  IF (J.EQ.4) KT = 229
                  IF (J.EQ.5) KT = 270
                  IF (J.EQ.6) KT = 279
               ENDIF
               IF(NPGSMX(I,J).EQ.0) GOTO 5
               DO K = 1,NPGSMX(I,J)
                  KT = KT + 1
                  MK = KT
                  IF (I.EQ.1 .AND. J.EQ.1 .AND. K.EQ.10) MK = 75
                  IF (I.EQ.1 .AND. J.EQ.3 .AND. K.GE.29) MK = -1
                  IF (I.EQ.1 .AND. J.EQ.4 .AND. K.GE.12) MK = -1
                  IF (I.EQ.1 .AND. J.EQ.5 .AND. K.GE. 5) MK = -1
                  IF (I.EQ.3 .AND. J.EQ.5 .AND. K.GE. 4) MK = -1
                  IF (I.EQ.3 .AND. J.EQ.6 .AND. K.GE. 3) MK = -1
                  IF (I.EQ.3 .AND. J.EQ.8 .AND. K.GE. 1) MK = -1
                  IF (I.EQ.3 .AND. J.EQ.9 .AND. K.GE. 1) MK = -1
                  MODECR(I,J,K) = MK
               ENDDO
    5          CONTINUE
            ENDDO
         ENDDO
C Define Printer
         IF(N6.LE.0) N6 = L6
C Define FIRST flag
         FIRST = .TRUE.
C Set the default values
         IF(IFLPRT.LT.0 .OR. IFLPRT.GT.5) IFLPRT = LFLPRT
         NPTYPE = LPTYPE
         NGROUP = LGROUP
         NSET   = LNSET
         NFL    = IFL
         LO     = LORD
         TMAS   = TM
         QCDL4  = PQCDL4(NPTYPE,NGROUP,NSET)
         QCDL5  = PQCDL5(NPTYPE,NGROUP,NSET)
         XMIN   = PXMIN(NPTYPE,NGROUP,NSET)
         XMAX   = PXMAX(NPTYPE,NGROUP,NSET)
         Q2MIN  = PQ2MIN(NPTYPE,NGROUP,NSET)
         Q2MAX  = PQ2MAX(NPTYPE,NGROUP,NSET)
C.
C Save the PDFLIB parameters
         IFLPRTP = IFLPRT
         NPTYPP  = NPTYPE
         NGROPP  = NGROUP
         NSETP   = NSET
         NFLP    = NFL
         LOP     = LO
         TMASP   = TMAS
         QCDL4P  = QCDL4
         QCDL5P  = QCDL5
         XMINP   = XMIN
         XMAXP   = XMAX
         Q2MINP  = Q2MIN
         Q2MAXP  = Q2MAX
C. endif from ISTART
      ENDIF
C.
C Validate call to PDFSET
      IFLSET = 1
C.
C Default :  NEW version 4.0 of PDFLIB
      NEWVER = .TRUE.
C.
C Now get user defined values
      NFLU   = 0
      LOU    = 0
      TMASU  = ZEROD
      QCDL4U = ZEROD
      QCDL5U = ZEROD
      XMINU  = ZEROD
      XMAXU  = ZEROD
      Q2MINU = ZEROD
      Q2MAXU = ZEROD
      KPTYP = 0
      KGROP = 0
      KNSET = 0
      KMODE = 0
      NV = 0
      DO 10 I = 1,20
      STRING = PARM(I)
      CALL CLTOU(STRING)
      NSCHECK = 0
      DO NS = 1,NSDIM
        IF(STRING(1:2).EQ.AV(NS)) NSCHECK = NS
      ENDDO
      IF(    STRING(1:5).NE.'INIT0' .AND. STRING(1:4).NE.'MODE'
     + .AND. STRING(1:6).NE.'NPTYPE'.AND. STRING(1:6).NE.'NGROUP'
     + .AND. STRING(1:4).NE.'NSET'  .AND. STRING(1:3).NE.'NFL'
     + .AND. STRING(1:2).NE.'LO'    .AND. STRING(1:4).NE.'TMAS'
     + .AND. STRING(1:5).NE.'QCDL4' .AND. STRING(1:5).NE.'QCDL5'
     + .AND. STRING(1:4).NE.'XMIN'  .AND. STRING(1:4).NE.'XMAX'
     + .AND. STRING(1:5).NE.'Q2MIN' .AND. STRING(1:5).NE.'Q2MAX'
C.
     + .AND. STRING(1:2).NE.'NU'    .AND. STRING(1:2).NE.'PR'
     + .AND. STRING(1:2).NE.'PI'    .AND. STRING(1:2).NE.'PH'
     + .AND. STRING(1:3).NE.'PRI'   .AND. NSCHECK    .EQ.0
     + .AND. STRING(1:2).NE.'OW'    .AND. STRING(1:2).NE.'SM'
     + .AND. STRING(1:4).NE.'GRVP'  .AND. STRING(1:5).NE.'GRV-P'
     + .AND. STRING(1:5).NE.'GRVPI' .AND. STRING(1:6).NE.'GRV-PI'
     + .AND. STRING(1:4).NE.'ABFK'  .AND. STRING(1:3).NE.'DOG'
     + .AND. STRING(1:4).NE.'DO-G'  .AND. STRING(1:4).NE.'DOPH'
     + .AND. STRING(1:5).NE.'DO-PH' .AND. STRING(1:2).NE.'DG'
     + .AND. STRING(1:2).NE.'LA'    .AND. STRING(1:2).NE.'GS'
     + .AND. STRING(1:3).NE.'GRV'   .AND. STRING(1:3).NE.'GRS'
     + .AND. STRING(1:4).NE.'WHIT'  .AND. STRING(1:3).NE.'SAS'
     + .AND. STRING(1:2).NE.'AC'    .AND. STRING(1:5).NE.'FIRST')GOTO 20
      NV = NV + 1
      VALUE  = VAL(I)
C.
      IF(STRING(1:5).EQ.'FIRST' .OR. STRING(1:4).EQ.'MODE'   .OR.
     +   STRING(1:6).EQ.'NPTYPE'.OR. STRING(1:6).EQ.'NGROUP' .OR.
     +   STRING(1:4).EQ.'NSET'  .OR. STRING(1:3).EQ.'NFL'    .OR.
     +   STRING(1:2).EQ.'LO'    .OR.
     +   STRING(1:3).EQ.'PRI'   .OR. NSCHECK    .GT.0        .OR.
     +   STRING(1:2).EQ.'OW'    .OR. STRING(1:2).EQ.'SM'     .OR.
     +   STRING(1:4).EQ.'GRVP'  .OR. STRING(1:5).EQ.'GRV-P'  .OR.
     +   STRING(1:5).EQ.'GRVPI' .OR. STRING(1:6).EQ.'GRV-PI' .OR.
     +   STRING(1:4).EQ.'ABFK'  .OR. STRING(1:3).EQ.'DOG'    .OR.
     +   STRING(1:4).EQ.'DO-G'  .OR. STRING(1:4).EQ.'DOPH'   .OR.
     +   STRING(1:5).EQ.'DO-PH' .OR. STRING(1:2).EQ.'DG'     .OR.
     +   STRING(1:2).EQ.'LA'    .OR. STRING(1:2).EQ.'GS'     .OR.
     +   STRING(1:3).EQ.'GRV'   .OR. STRING(1:3).EQ.'GRS'    .OR.
     +   STRING(1:4).EQ.'WHIT'  .OR. STRING(1:3).EQ.'SAS'    .OR.
     +   STRING(1:2).EQ.'AC'                               ) THEN
         TEMP = VALUE
         IF(ABS(VALUE).LT.1.E-2) TEMP = IVALUE
      ENDIF
C Get user's setting for FIRST (option: silent mode => FIRST = .FALSE.)
      IF(STRING(1:5).EQ.'FIRST')   THEN
        IF(TEMP.EQ.0) FIRST = .FALSE.
        IF(TEMP.EQ.1) FIRST = .TRUE.
      ENDIF
C.
C. Return, because PDFSET was called ONLY for initialisation
      IF(STRING(1:5).EQ.'INIT0') THEN
         IFLSET = 0
C Print status of PDFLIB
         IF(FIRST) THEN
            IF(ISTART.EQ.0) CALL PDFVERS
         ENDIF
         ISTART = 1
C. Reset PDFSET parameters
         DO K = 1,NCHDIM
            PARM(K) = '  '
         ENDDO
         RETURN
      ENDIF
C.
C Check if PDFLIB format of new version 4.0 or of versions 3.0 or less
      IF(STRING(1:4).EQ.'MODE') THEN
        NEWVER = .FALSE.
        MODE = TEMP
        KMODE = 1
      ENDIF
C Now NEW version 4.0 or bigger
C Get particle type
      IF(STRING(1:6).EQ.'NPTYPE')   THEN
        NPTYPE = TEMP
        KPTYP  = 1
      ENDIF
C Get Author Group Number
      IF(STRING(1:6).EQ.'NGROUP') THEN
        NGROUP = TEMP
        KGROP  = 1
      ENDIF
C Get PDF Set within Author Group
      IF(STRING(1:4).EQ.'NSET' ) THEN
         NSET  = TEMP
         KNSET = 1
      ENDIF
C Particle Type = Nucleons
      IF(STRING(1:2).EQ.'NU' .OR. STRING(1:2).EQ.'PR'
     +                       .OR. NSCHECK.GT.0)        THEN
        NPTYPE = 1
        KPTYP  = 1
      ENDIF
      IF(STRING(1:3).EQ.'PRI') THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 1
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.GE.1 .AND. NSCHECK.LE.5) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 1
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.EQ.6) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 2
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.GE.7 .AND. NSCHECK.LE.9) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 3
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.GE.10 .AND. NSCHECK.LE.11) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 4
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.EQ.12) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 5
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.EQ.13) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 6
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(NSCHECK.EQ.14) THEN
        NPTYPE = 1
        KPTYP  = 1
        NGROUP = 7
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
C Particle Type = Pions
      IF(STRING(1:2).EQ.'PI') THEN
        NPTYPE = 2
        KPTYP  = 1
      ENDIF
      IF(STRING(1:2).EQ.'OW'   .OR.
     +   STRING(1:3).EQ.'OWP'  .OR. STRING(1:4).EQ.'OW-P'   .OR.
     +   STRING(1:4).EQ.'OWPI' .OR. STRING(1:5).EQ.'OW-PI') THEN
        NPTYPE = 2
        KPTYP  = 1
        NGROUP = 1
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'SM'     .OR.
     +   STRING(1:5).EQ.'SMRSP'  .OR. STRING(1:6).EQ.'SMRS-P'   .OR.
     +   STRING(1:6).EQ.'SMRSPI' .OR. STRING(1:7).EQ.'SMRS-PI') THEN
        NPTYPE = 2
        KPTYP  = 1
        NGROUP = 3
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:4).EQ.'GRVP'  .OR. STRING(1:5).EQ.'GRV-P'   .OR.
     +   STRING(1:5).EQ.'GRVPI' .OR. STRING(1:6).EQ.'GRV-PI') THEN
        NPTYPE = 2
        KPTYP  = 1
        NGROUP = 5
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:4).EQ.'ABFK'    .OR.
     +   STRING(1:6).EQ.'ABFKWP'  .OR. STRING(1:7).EQ.'ABFKW-P'   .OR.
     +   STRING(1:7).EQ.'ABFKWPI' .OR. STRING(1:8).EQ.'ABFKW-PI') THEN
        NPTYPE = 2
        KPTYP  = 1
        NGROUP = 6
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
C Particle Type = Photons
      IF(STRING(1:2).EQ.'PH') THEN
        NPTYPE = 3
        KPTYP  = 1
      ENDIF
      IF(STRING(1:3).EQ.'DOG'  .OR. STRING(1:4).EQ.'DO-G'   .OR.
     +   STRING(1:4).EQ.'DOPH' .OR. STRING(1:5).EQ.'DO-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 1
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'DG'   .OR.
     +   STRING(1:3).EQ.'DGG'  .OR. STRING(1:4).EQ.'DG-G'   .OR.
     +   STRING(1:4).EQ.'DGPH' .OR. STRING(1:5).EQ.'DG-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 2
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'LA'    .OR.
     +   STRING(1:4).EQ.'LACG'  .OR. STRING(1:5).EQ.'LAC-G'   .OR.
     +   STRING(1:5).EQ.'LACPH' .OR. STRING(1:6).EQ.'LAC-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 3
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'GS'   .OR.
     +   STRING(1:3).EQ.'GSG'  .OR. STRING(1:4).EQ.'GS-G'   .OR.
     +   STRING(1:4).EQ.'GSPH' .OR. STRING(1:5).EQ.'GS-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 4
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:4).EQ.'GRVG'  .OR. STRING(1:5).EQ.'GRV-G'   .OR.
     +   STRING(1:5).EQ.'GRVPH' .OR. STRING(1:6).EQ.'GRV-PH'  .OR.
     +   STRING(1:5).EQ.'GRSPH' .OR. STRING(1:6).EQ.'GRS-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 5
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'AC'      .OR.
     +   STRING(1:6).EQ.'ACFGPG'  .OR. STRING(1:7).EQ.'ACFGP-G'   .OR.
     +   STRING(1:7).EQ.'ACFGPPH' .OR. STRING(1:8).EQ.'ACFGP-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 6
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'WI'      .OR.
     +   STRING(1:5).EQ.'WHITG'  .OR. STRING(1:6).EQ.'WHIT-G'   .OR.
     +   STRING(1:6).EQ.'WHITPH' .OR. STRING(1:7).EQ.'WHIT-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 8
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
      IF(STRING(1:2).EQ.'SA'      .OR.
     +   STRING(1:4).EQ.'SASG'  .OR. STRING(1:5).EQ.'SAS-G'   .OR.
     +   STRING(1:5).EQ.'SASPH' .OR. STRING(1:6).EQ.'SAS-PH') THEN
        NPTYPE = 3
        KPTYP  = 1
        NGROUP = 9
        KGROP  = 1
        NSET   = TEMP
        KNSET  = 1
      ENDIF
C... set defaults for private structure functions set (Nset = 1)
      IF(NPTYPE.EQ.1 .AND. NGROUP.EQ.1 .AND. NSET.EQ.1) THEN
         NFL   = IFL
         LO    = LORD
         TMAS  = TM
         QCDL4 = ZEROD
         QCDL5 = ZEROMIN
         XMIN  = ZEROD
         XMAX  = 0.9999D0
         Q2MIN = ZEROD
         Q2MAX = 1.000D8
      ENDIF
      IF(STRING(1:3).EQ.'NFL'  ) NFLU   = TEMP
      IF(STRING(1:2).EQ.'LO'   ) LOU    = TEMP
      IF(STRING(1:4).EQ.'TMAS' ) TMASU  = VALUE
C.
      IF(STRING(1:5).EQ.'QCDL4') QCDL4U = VALUE
      IF(STRING(1:5).EQ.'QCDL5') QCDL5U = VALUE
C.
      IF(STRING(1:4).EQ.'XMIN')  XMINU  = VALUE
      IF(STRING(1:4).EQ.'XMAX')  XMAXU  = VALUE
      IF(STRING(1:5).EQ.'Q2MIN') Q2MINU = VALUE
      IF(STRING(1:5).EQ.'Q2MAX') Q2MAXU = VALUE
   10 CONTINUE
   20 CONTINUE
C Print status of PDFLIB
         IF(FIRST) THEN
            IF(ISTART.EQ.0) CALL PDFVERS
         ENDIF
         ISTART = 1
C Print values
      IF(IFLPRT.GE.1) THEN
         IF(NV.GT.0) THEN
            IF(NEWVER) THEN
               WRITE(N6,2000) (PARM(I)(1:6),I=1,NV)
            ELSE
               WRITE(N6,2000) (PARM(I)(1:4),I=1,NV)
            ENDIF
            WRITE(N6,2001) (VAL(I),I=1,NV)
         ELSE
            WRITE(N6,2002)
         ENDIF
      ELSE
         IF(NV.EQ.0) WRITE(N6,2002)
      ENDIF
 2000 FORMAT(/,'  Parm = ',10(1X,A8))
 2001 FORMAT(1H ,' Val  = ',10(1X,F8.4))
 2002 FORMAT(/,'  PDFSET : Warning!!  NO Parameters read in,   defaults
     + (Nucleon PDFs) used !!')
C.
C. check on conflicting parameters in case of mixed versions
      IF((KPTYP*KMODE).EQ.1 .OR. (KGROP*KMODE).EQ.1 .OR.
     +   (KNSET*KMODE).EQ.1) THEN
           WRITE(N6,*)
     +   ' PDFLIB  : Warning!!   Mixed parameter settings for OLD and NE
     +W version of PDFLIB format'
           WRITE(N6,*)
     +  ' Warning : ALL parameters set to NEW default (Nucleon PDFs) !!'
           NEWVER = .TRUE.
           NPTYPE = LPTYPE
           NGROUP = LGROUP
           NSET   = LNSET
      ENDIF
C.
#include "pdf/w50511c.inc"
C.
C User wants old version (3.0 or less) of PDFLIB format
      IF(.NOT.NEWVER) THEN
        IF (MODE .GE.0 .AND. MODE .LE.MODEMX) THEN
           IF(MODE.EQ.0) THEN
             NPTYPE = LPTYPO
             NGROUP = LGROPO
             NSET   = LNSETO
           ELSE
             NPTYPE = NPTYCR(MODE)
             NGROUP = NGROCR(MODE)
             NSET   = NSETCR(MODE)
C Check on validity of parameter values
             IF (NPTYPE.LT.0 .OR. NGROUP.LT.0 .OR. NSET.LT.0) THEN
                IF(FIRST) THEN
                WRITE(N6,*) ' PDFLIB  : MODE value INCORRECT, MODE = ',
     +          MODE
               WRITE(N6,*) ' Warning : MODE value set to OLD default !!'
                ENDIF
                NPTYPE = LPTYPO
                NGROUP = LGROPO
                NSET   = LNSETO
                GOTO 30
             ENDIF
          ENDIF
        ELSE
          IF(FIRST) THEN
          WRITE(N6,*) ' PDFLIB  : MODE value INCORRECT, MODE = ',MODE
          WRITE(N6,*) ' Warning : MODE value set to OLD default !!'
          ENDIF
          NPTYPE = LPTYPO
          NGROUP = LGROPO
          NSET   = LNSETO
        ENDIF
C User wants new version (4.0 or bigger) of PDFLIB format
      ELSE
C Check on validity of parameter values
        IF (NSET.EQ.0) THEN
           NPTYPE = LPTYPE
           NGROUP = LGROUP
           NSET   = LNSET
           IF(FIRST) THEN
           WRITE(N6,*)
     +      ' PDFLIB  : DEFAULT parameter set (Nucleon PDFs) has been
     + chosen !'
           WRITE(N6,*)
     +      ' PDFLIB  : Nptype = ',NPTYPE,'   Ngroup = ',NGROUP,
     +      '   Nset = ',NSET
           ENDIF
        ENDIF
        IF (NPTYPE.LT.0 .OR. NPTYPE.EQ.0 .OR. NPTYPE.GT.NPTYMX) THEN
           IF(FIRST) THEN
           WRITE(N6,*)
     +   ' PDFLIB  : Particle type number INCORRECT, Nptype = ',
     +     NPTYPE,'   Ngroup = ',NGROUP,'   Nset = ',NSET
           WRITE(N6,*)
     +   ' Warning : ALL parameters set to default (Nucleon PDFs) !!'
           ENDIF
           NPTYPE = LPTYPE
           NGROUP = LGROUP
           NSET   = LNSET
        ENDIF
        IF (NGROUP.LT.0 .OR. NGROUP.EQ.0 .OR. NGROUP.GT.NGRMAX) THEN
           IF(FIRST) THEN
           WRITE(N6,*)
     +   ' PDFLIB  : Group number for PDFs INCORRECT, Nptype = ',
     +     NPTYPE,'   Ngroup = ',NGROUP,'   Nset = ',NSET
           WRITE(N6,*)
     +   ' Warning : ALL parameters set to default (Nucleon PDFs) !!'
           ENDIF
           NPTYPE = LPTYPE
           NGROUP = LGROUP
           NSET   = LNSET
        ENDIF
        IF (NSET.LT.0 .OR. NSET.EQ.0 .OR. NSET.GT.NPGSMX(NPTYPE,NGROUP))
     +  THEN
           IF(FIRST) THEN
           WRITE(N6,*)
     +   ' PDFLIB  : NSET value for PDFs INCORRECT, Nptype = ',
     +     NPTYPE,'   Ngroup = ',NGROUP,'   Nset = ',NSET
           WRITE(N6,*)
     +   ' Warning : ALL parameters set to default (Nucleon PDFs) !!'
           ENDIF
           NPTYPE = LPTYPE
           NGROUP = LGROUP
           NSET   = LNSET
        ENDIF
      ENDIF
   30 CONTINUE
C.
C Check on retracted PDF sets
        IF (NSETFL(NPTYPE,NGROUP,NSET).EQ.0) THEN
           IF(FIRST) THEN
           WRITE(N6,*)
     +   ' PDFLIB  : Your selected PDF set has been retracted by the aut
     +hors'
           WRITE(N6,*)
     +   ' Warning : ALL parameters set to default (Nucleon PDFs) !!'
           ENDIF
           NPTYPE = LPTYPE
           NGROUP = LGROUP
           NSET   = LNSET
        ENDIF
C.
C Number of flavours in alpha(s) calculation
      IF (ABS(NFL).LT.3 .OR. ABS(NFL).GT.6) THEN
         IF(FIRST) THEN
         IF(NPTYPE.EQ.1 .AND. NGROUP.EQ.1 .AND. NSET.EQ.1) WRITE(N6,*)
     +   '   You are calling the DUMMY private SF routine !!'
         WRITE(N6,*) ' PDFLIB  : NFL value outside range, NFL = ',NFL
         WRITE(N6,*) ' Warning : NFL value set to default !!'
         ENDIF
         NFL = IFL
      ENDIF
      IF (NFLU .GT. 0) THEN
         NFL = NFLU
         IF(FIRST) THEN
         IF(NFLU.NE.NFLP) THEN
         WRITE(N6,*) ' PDFLIB  : NFL value'
         WRITE(N6,*) ' Warning : NON standard settings, NFL value = ',
     +   NFL,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Order of alpha(s) calculation
      LO=LORD
      IF ((NPTYPE.EQ.1.AND.NGROUP.EQ.1.AND.NSET.GE. 2) .OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.2.AND.NSET.LE. 5) .OR.
     +   (NPTYPE.EQ.1.AND.NGROUP.EQ.4.AND.NSET.GE.11.AND.NSET.LE.12).OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.4.AND.NSET.EQ.17) .OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.4.AND.NSET.EQ.23) .OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.4.AND.NSET.EQ.29) .OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.4.AND.NSET.EQ.32) .OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.5.AND.NSET.EQ. 2) .OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.5.AND.NSET.GE. 4.AND.NSET.LE.5).OR.
     +    (NPTYPE.EQ.1.AND.NGROUP.EQ.5.AND.NSET.GE. 8.AND.NSET.LE.9).OR.
     +    (NPTYPE.EQ.2.AND.NGROUP.EQ.1.AND.NSET.LE. 2) .OR.
     +    (NPTYPE.EQ.2.AND.NGROUP.EQ.5.AND.NSET.EQ. 2) .OR.
     +    (NPTYPE.EQ.3.AND.NGROUP.EQ.1.AND.NSET.EQ. 1) .OR.
     +    (NPTYPE.EQ.3.AND.NGROUP.EQ.2.AND.NSET.LE. 4) .OR.
     +    (NPTYPE.EQ.3.AND.NGROUP.EQ.3.AND.NSET.LE. 3) .OR.
     +    (NPTYPE.EQ.3.AND.NGROUP.EQ.4.AND.NSET.GE. 2.AND.NSET.LE.3).OR.
     +    (NPTYPE.EQ.3.AND.NGROUP.EQ.4.AND.NSET.EQ. 5) .OR.
     +    (NPTYPE.EQ.3.AND.NGROUP.EQ.5.AND.NSET.GE. 3.AND.NSET.LE.4).OR.
     +    (NPTYPE.EQ.3.AND.(NGROUP.GE.8.AND.NGROUP.LE.9)))
     +     LO=1
      IF (LOU .GT. 0) THEN
         LO = LOU
         IF(FIRST) THEN
         IF(LOU.NE.LOP) THEN
         WRITE(N6,*) ' PDFLIB  : LO  value'
         WRITE(N6,*) ' Warning : NON standard settings, LO value = ',
     +   LO,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Top Quark mass
      IF(TMAS.LE.ZEROD) THEN
        IF(FIRST) THEN
        WRITE(N6,*)' PDFLIB  : TMAS value zero or negativ, Tmas = ',TMAS
        WRITE(N6,*)' Warning : TMAS value set to default !!'
        ENDIF
        TMAS = TM
      ENDIF
      IF (TMASU .GT. ZEROMIN) THEN
         TMAS = TMASU
         IF(FIRST) THEN
         IF(TMASU.NE.TMASP) THEN
         WRITE(N6,*)' PDFLIB  : TMAS value'
         WRITE(N6,*) ' Warning : NON standard settings, TMAS value = ',
     +   TMAS,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C. this lambda (NFL=4) is given by the authors and
C.                     NOT used in the calculation
      QCDL4 = PQCDL4(NPTYPE,NGROUP,NSET)
      IF(NPTYPE.EQ.1 .AND. NGROUP.EQ.1 .AND. NSET.EQ.1) THEN
         IF (QCDL5 .LT. ZEROMIN) THEN
            IF(FIRST) THEN
            WRITE (N6,*) ' PDFLIB  : QCD-Lambda5 is not defined !! '
            WRITE (N6,*) ' Warning : QCD-Lambda5 set to QCDL5 = 0.0D0'
            ENDIF
            QCDL4 = ZEROD
            QCDL5 = ZEROD
C. Reset PDFSET parameters
            DO K = 1,NCHDIM
               PARM(K) = '  '
            ENDDO
            RETURN
         ENDIF
      ENDIF
      IF (QCDL4U .GT. ZEROMIN) THEN
         QCDL4 = QCDL4U
         IF(FIRST) THEN
         IF(QCDL4U.NE.QCDL4P) THEN
         WRITE (N6,*) ' PDFLIB  : QCD-Lambda4 value'
         WRITE(N6,*) ' Warning : NON standard settings, QCDL4 value = ',
     +   QCDL4,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C. this is the lambda (NFL=5) used in the calculation
      QCDL5 = PQCDL5(NPTYPE,NGROUP,NSET)
      IF (QCDL5U .GT. ZEROMIN) THEN
         QCDL5 = QCDL5U
         IF(FIRST) THEN
         IF(QCDL5U.NE.QCDL5P) THEN
         WRITE (N6,*) ' PDFLIB  : QCD-Lambda5 value'
         WRITE(N6,*) ' Warning : NON standard settings, QCDL5 value = ',
     +   QCDL5,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Minimal X values
      XMIN  = PXMIN(NPTYPE,NGROUP,NSET)
      IF (XMINU .GT. ZEROMIN) THEN
         XMIN = XMINU
         IF(FIRST) THEN
         IF(XMINU.NE.XMINP) THEN
         WRITE (N6,*) ' PDFLIB  : XMIN value'
         WRITE(N6,*) ' Warning : NON standard settings, XMIN value = ',
     +   XMIN,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Maximal X values
      XMAX  = PXMAX(NPTYPE,NGROUP,NSET)
      IF (XMAXU .GT. ZEROMIN) THEN
         XMAX = XMAXU
         IF(FIRST) THEN
         IF(XMAXU.NE.XMAXP) THEN
         WRITE (N6,*) ' PDFLIB  : XMAX value'
         WRITE(N6,*) ' Warning : NON standard settings, XMAX value = ',
     +   XMAX,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Minimal Q**2 values
      Q2MIN = PQ2MIN(NPTYPE,NGROUP,NSET)
      IF (Q2MINU .GT. ZEROMIN) THEN
         Q2MIN = Q2MINU
         IF(FIRST) THEN
         IF(Q2MINU.NE.Q2MINP) THEN
         WRITE (N6,*) ' PDFLIB  : Q2MIN value'
         WRITE(N6,*) ' Warning : NON standard settings, Q2MIN value = ',
     +   Q2MIN,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Maximal Q**2 values
      Q2MAX = PQ2MAX(NPTYPE,NGROUP,NSET)
      IF (Q2MAXU .GT. ZEROMIN) THEN
         Q2MAX = Q2MAXU
         IF(FIRST) THEN
         IF(Q2MAXU.NE.Q2MAXP) THEN
         WRITE (N6,*) ' PDFLIB  : Q2MAX value'
         WRITE(N6,*) ' Warning : NON standard settings, Q2MAX value = ',
     +   Q2MAX,'  set by user !!'
         ENDIF
         ENDIF
      ENDIF
C.
C Save the PDFLIB parameters
         IFLPRTP = IFLPRT
         NPTYPP  = NPTYPE
         NGROPP  = NGROUP
         NSETP   = NSET
         NFLP    = NFL
         LOP     = LO
         TMASP   = TMAS
         QCDL4P  = QCDL4
         QCDL5P  = QCDL5
         XMINP   = XMIN
         XMAXP   = XMAX
         Q2MINP  = Q2MIN
         Q2MAXP  = Q2MAX
C.
C... Reset variables to zero for error counting
           PDFWGT = ONED
           WXMIN = ZEROD
           WXMAX = ZEROD
           WQ2MIN = ZEROD
           WQ2MAX = ZEROD
           WTXMIN = ZEROD
           WTXMAX = ZEROD
           WTQ2MIN = ZEROD
           WTQ2MAX = ZEROD
C.
C print COMMON block values
      IF(IFLPRT.GE.2) THEN
         WRITE(N6,3000) NPTYPE,NGROUP,NSET,SFNAME(NPTYPE,NGROUP,NSET),
     +                  MODECR(NPTYPE,NGROUP,NSET)
         WRITE(N6,3001) NFL,LO,TMAS
         WRITE(N6,3002) QCDL4,QCDL5
         WRITE(N6,3003) XMIN,XMAX,Q2MIN,Q2MAX
      ENDIF
 3000 FORMAT(/,'  Nptype = ',I1,'  Ngroup = ',I1,'  Nset = ',I2,
     +         '  Name  = "',A8,'"  CrMode = ',I3)
 3001 FORMAT(1H ,' Nfl    = ',I2,',  LO = ',I1,
     +           ',  Tmas = ',F7.2,' GeV/c**2')
 3002 FORMAT(1H ,' QCDL4  = ',F7.4,' GeV,  QCDL5 = ',F7.4,' GeV')
 3003 FORMAT(1H ,' Xmin   = ',1X,E8.2,',  Xmax = ',E11.5,
     +          ',  Q2min = ',F7.3,' (GeV/c)**2,  Q2max = ',E8.2,
     +          ' (GeV/c)**2')
C.
C. Reset PDFSET parameters
      DO K = 1,NCHDIM
         PARM(K) = '  '
      ENDDO
C.
      RETURN
      END

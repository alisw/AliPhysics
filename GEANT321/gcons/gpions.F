*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GPIONS
C.
C.    *******************************************************************
C.    *                                                                 *
C.    *  Define a subset of the 'stable' most common elements           *
C.    *    in the Nature                                                *
C.    *                                                                 *
C.    *    ==>Called by : <USER>, UGINIT                                *
C.    *       Author    B.Grosdidier (Strasbourg) *****                 *
C.    *                                                                 *
C.    *******************************************************************
C.
      DIMENSION UBUF(1)
      DATA UBUF/0./
*
*     The Geant particle identification begin to IPART = 61
      CALL GSPART( 61,'LI6  ',8,   5.60305,  3., 1000., UBUF,0)
      CALL GSPART( 62,'LI7  ',8,   6.53536,  3., 1000., UBUF,0)
      CALL GSPART( 63,'BE7  ',8,   6.53622,  4., 1000., UBUF,0)
      CALL GSPART( 64,'BE9  ',8,   8.39479,  4., 1000., UBUF,0)
      CALL GSPART( 65,'B10  ',8,   9.32699,  5., 1000., UBUF,0)
      CALL GSPART( 66,'B11  ',8,  10.25510,  5., 1000., UBUF,0)
      CALL GSPART( 67,'C12  ',8,  11.17793,  6., 1000., UBUF,0)
      CALL GSPART( 68,'N14  ',8,  13.04378,  7., 1000., UBUF,0)
      CALL GSPART( 69,'O16  ',8,  14.89917,  8., 1000., UBUF,0)
      CALL GSPART( 70,'F19  ',8,  17.69690,  9., 1000., UBUF,0)
      CALL GSPART( 71,'NE20 ',8,  18.62284, 10., 1000., UBUF,0)
      CALL GSPART( 72,'NA23 ',8,  21.41483, 11., 1000., UBUF,0)
      CALL GSPART( 73,'MG24 ',8,  22.34193, 12., 1000., UBUF,0)
      CALL GSPART( 74,'AL27 ',8,  25.13314, 13., 1000., UBUF,0)
      CALL GSPART( 75,'SI28 ',8,  26.06034, 14., 1000., UBUF,0)
      CALL GSPART( 76,'P31  ',8,  28.85188, 15., 1000., UBUF,0)
      CALL GSPART( 77,'S32  ',8,  29.78180, 16., 1000., UBUF,0)
      CALL GSPART( 78,'CL35 ',8,  32.57328, 17., 1000., UBUF,0)
      CALL GSPART( 79,'AR36 ',8,  33.50356, 18., 1000., UBUF,0)
      CALL GSPART( 80,'K39  ',8,  36.29447, 19., 1000., UBUF,0)
      CALL GSPART( 81,'CA40 ',8,  37.22492, 20., 1000., UBUF,0)
      CALL GSPART( 82,'SC45 ',8,  41.87617, 21., 1000., UBUF,0)
      CALL GSPART( 83,'TI48 ',8,  44.66324, 22., 1000., UBUF,0)
      CALL GSPART( 84,'V51  ',8,  47.45401, 23., 1000., UBUF,0)
      CALL GSPART( 85,'CR52 ',8,  48.38228, 24., 1000., UBUF,0)
      CALL GSPART( 86,'MN55 ',8,  51.17447, 25., 1000., UBUF,0)
      CALL GSPART( 87,'FE56 ',8,  52.10307, 26., 1000., UBUF,0)
      CALL GSPART( 88,'CO59 ',8,  54.89593, 27., 1000., UBUF,0)
      CALL GSPART( 89,'NI58 ',8,  53.96644, 28., 1000., UBUF,0)
      CALL GSPART( 90,'CU63 ',8,  58.61856, 29., 1000., UBUF,0)
      CALL GSPART( 91,'ZN64 ',8,  59.54963, 30., 1000., UBUF,0)
      CALL GSPART( 92,'GE74 ',8,  68.85715, 32., 1000., UBUF,0)
      CALL GSPART( 93,'SE80 ',8,  74.44178, 34., 1000., UBUF,0)
      CALL GSPART( 94,'KR84 ',8,  78.16309, 36., 1000., UBUF,0)
      CALL GSPART( 95,'SR88 ',8,  81.88358, 38., 1000., UBUF,0)
      CALL GSPART( 96,'ZR90 ',8,  83.74571, 40., 1000., UBUF,0)
      CALL GSPART( 97,'MO98 ',8,  91.19832, 42., 1000., UBUF,0)
      CALL GSPART( 98,'PD106',8,  98.64997, 46., 1000., UBUF,0)
      CALL GSPART( 99,'CD114',8, 106.10997, 48., 1000., UBUF,0)
      CALL GSPART(100,'SN120',8, 111.68821, 50., 1000., UBUF,0)
      CALL GSPART(101,'XE132',8, 122.86796, 54., 1000., UBUF,0)
      CALL GSPART(102,'BA138',8, 128.45793, 56., 1000., UBUF,0)
      CALL GSPART(103,'CE140',8, 130.32111, 58., 1000., UBUF,0)
      CALL GSPART(104,'SM152',8, 141.51236, 62., 1000., UBUF,0)
      CALL GSPART(105,'DY164',8, 152.69909, 66., 1000., UBUF,0)
      CALL GSPART(106,'YB174',8, 162.02245, 70., 1000., UBUF,0)
      CALL GSPART(107,'W184 ',8, 171.34924, 74., 1000., UBUF,0)
      CALL GSPART(108,'PT194',8, 180.67513, 78., 1000., UBUF,0)
      CALL GSPART(109,'AU197',8, 183.47324, 79., 1000., UBUF,0)
      CALL GSPART(110,'HG202',8, 188.13451, 80., 1000., UBUF,0)
      CALL GSPART(111,'PB208',8, 193.72907, 82., 1000., UBUF,0)
      CALL GSPART(112,'U238 ',8, 221.74295, 92., 1000., UBUF,0)
*
      END

*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:15:28  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:46:56  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_TZC1_INC
#define CERNLIB_ZEBRA_TZC1_INC
*
*
* tzc1.inc
*
                     PARAMETER  (LGL=240)
                     CHARACTER   LINE*(LGL), LHEAD*80, BLANK*(LGL)
      COMMON /TZC1/  LUNTQ,  NREADY, NFATAL, IGNORE
     +,              IFLLOG, IFLPRI, IFLUSE, IFLKIL, IFLANY
     +,              IFLEXA, IFLSIZ, IFLTZU, IVALUS, JFMTC1, JFMTC2
     +,              JCOLA,  JCOLE,  NCHPW,  IFLHC,  IFLHW,  IFLDBL
     +,              IFLPRE, NFAULT, LPUTA,  LPUTX,  LPUTE
     +,              NCHORG, NHEAD,  LINE,  LHEAD,  BLANK

#endif

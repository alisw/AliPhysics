*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:53  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_MZCT_INC
#define CERNLIB_ZEBRA_MZCT_INC
*
*                 Memory Table parameters
*
* mzct.inc
*
      COMMON /MZCT/  MQDVGA,MQDVWI,JQSTMV,JQDVM1,JQDVM2,NQDVMV,IQFLIO
     +,              MQDVAC,NQNOOP,IQPART,NQFREE, IQTBIT,IQTVAL
     +,              IQTNMV,JQGAPM,JQGAPR,NQGAPN,NQGAP,IQGAP(5,4)
     +,              LQTA,LQTE, LQRTA,LQTC1,LQTC2,LQRTE
     +,                         LQMTA,LQMTB,LQMTE,LQMTLU,LQMTBR
     +,              LQMTC1,LQMTC2, NQFRTC,NQLIVE

#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_ZKRAKCC_INC
#define CERNLIB_ZEBRA_ZKRAKCC_INC
*
*
* zkrakcc.inc
*
#if defined(CERNLIB_EQUHOLCH)
      COMMON /ZKRAKC/IQHOLK(120), IQKRAK(80), IQCETK(122)
                     CHARACTER    CQHOLK*120
                     EQUIVALENCE (CQHOLK,IQHOLK)
#else
      COMMON /ZKRAKC/IQHOLK(120), IQKRAK(80), IQCETK(122)
      COMMON /ZCRAKC/CQHOLK
                     CHARACTER    CQHOLK*120
#endif

#endif

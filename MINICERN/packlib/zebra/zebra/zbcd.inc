*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_ZBCD_INC
#define CERNLIB_ZEBRA_ZBCD_INC
*
*                 BCD char. in CETA order with A1 format
*
* zbcd.inc
*
      COMMON /ZBCD/  IQNUM2(11),IQLETT(26),IQNUM(10),   IQPLUS,IQMINS
     +,              IQSTAR,IQSLAS,IQOPEN,IQCLOS,IQDOLL,IQEQU, IQBLAN
     +,              IQCOMA,IQDOT, IQNUMB,IQAPO, IQEXCL,IQCOLO,IQQUOT
     +,              IQUNDE,IQCLSQ,IQAND, IQAT,  IQQUES,IQOPSQ,IQGREA
     +,              IQLESS,IQREVE,IQCIRC,IQSEMI,IQPERC,  IQLOWL(26)
     +,              IQCROP,IQVERT,IQCRCL,IQNOT, IQGRAV,  IQILEG
     +,              NQHOL0,NQHOLL(95)

#endif

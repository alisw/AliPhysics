*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:55  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_HZC_INC
#define CERNLIB_ZEBRA_HZC_INC
*
*     HISTOGRAM MANIPULATION
*
* hzc.inc
*
#if !defined(CERNLIB_QCDE)
      COMMON /HZC/   NQHNH,NQHPAK(3),NQHGW,NQHNX,QHDX,QHXL,QHSTEP,NQHAND
     +,              LQHH,LQHE,LQHC,LQHL,LQHLB,LQHLC,LQHD,LQHX,LQHT,LQHW
     +,              IQHCUT,IQHXXX,IQHINT,NQHLAB,IQHFMT,IQHW(8)
     +,              QHSMA,QHBIG,QHTOT,NQHNX2
     +,              LQHCUT,NQHEXP,IQHOPT(7)
*  -------------   for defaults   ------------------------------
#endif

#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*
* fzudat49.inc
*
#if defined(CERNLIB_QMIBMD)
      CALL ADWRTKD (LUNX, LQ(LWR),NWMIN,1,ISW)
*    --------------  Here we are , MRo        ------------------
#endif
*    --------------  conversions for input   -------------------
*
* fzudat49.inc
*
#if defined(CERNLIB_QMIRTD)
      CALL ADWRTKD (LUNX, LQ(LWR),NWMIN,1,ISW)
*    --------------  conversions for output  -------------------
#endif
*    use default CALL IE3FOS for output single prec.
#ifndef CERNLIB_FQIE3FSC
#endif
*    use default CALL IE3FOD for output double prec.
#ifndef CERNLIB_FQIE3FDC
#endif
*    --------------  sequences for input   ---------------------

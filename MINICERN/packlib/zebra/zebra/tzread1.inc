*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/18 16:15:29  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:46:56  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_TZREAD1_INC
#define CERNLIB_ZEBRA_TZREAD1_INC
*
*
* tzread1.inc
*
#if defined(CERNLIB_QMIBM)
      READ (LUNTQ,NUM=NCH,END=78) LINE
      NCH    = MIN (NCH,LGL)
      NCHORG = LNBLNK (LINE(1:NCH))
#elif defined(CERNLIB_QMIBMD)
      READ (LUNTQ,NUM=NCH,END=78) LINE
      NCH    = MIN (NCH,LGL)
      NCHORG = LNBLNK (LINE(1:NCH))
#elif defined(CERNLIB_QMVAX)
      READ (LUNTQ,8000,END=78) NCH,LINE(1:MIN(NCH,LGL))
      NCH    = MIN (NCH,LGL)
      NCHORG = LNBLNK (LINE(1:NCH))
 8000 FORMAT (Q,A)
#elif 1
      READ (LUNTQ,'(A)',END=78) LINE
      NCHORG = LENOCC (LINE)
#endif

#endif

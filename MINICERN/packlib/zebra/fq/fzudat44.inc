*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*   pack 32-bit words to 60-bit
*
* fzudat44.inc
*
#if defined(CERNLIB_QMCDC)
      CALL PKCHAR (LQ(LWR),LQ(LWR-1),MAXREX,IPKVEC)
      LWR  = LWR - 1
#elif defined(CERNLIB_QMCRY)
      CALL PACK (LQ(LWR),32,LQ(LWR),MAXREX)
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      CALL PACK32 (LQ(LWR),LQ(LWR),MAXREX)
#elif defined(CERNLIB_QMIBMD)
      CALL PACK (LQ(LWR),LQ(LWR),MAXREX)
#elif defined(CERNLIB_QMIRTD)
      CALL PACK (LQ(LWR),LQ(LWR),MAXREX)
#elif 1
      CALL VXINVB (LQ(LWR),24)
#endif

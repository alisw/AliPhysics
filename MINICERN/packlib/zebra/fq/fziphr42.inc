*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:09  mclareni
* Zebra
*
*
*   unpack n-bit words to 32-bit
*
* fziphr42.inc
*
#if defined(CERNLIB_QMCDC)
   42 CALL UPKCH (LQ(LIN),LQ(L4STAI),NW4USE,IPKVEC)
#elif defined(CERNLIB_QMCRY)
   42 CALL UNPACK (LQ(LIN),32,LQ(L4STAI),NW4USE)
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
   42 CALL UPAK32 (LQ(LIN),LQ(L4STAI),NW4USE)
#elif defined(CERNLIB_QMIBMD)
   42 CALL UPACK (LQ(LIN),LQ(L4STAI),NW4USE)
#elif defined(CERNLIB_QMIRTD)
   42 CALL UPACK (LQ(LIN),LQ(L4STAI),NW4USE)
#elif 1
      NSWDON = 0
   42 CALL VXINVB (LQ(LIN+NSWDON),NW4USE-NSWDON)
      NSWDON = NW4USE
#endif

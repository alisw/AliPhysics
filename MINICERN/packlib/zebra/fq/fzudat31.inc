*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*   unpack 32-bit words to 60-bit
*
* fzudat31.inc
*
#if defined(CERNLIB_QMCDC)
      CALL UPKCH (LQ(LIN),LQ(L4STOX),MAXREX,IPKVEC)
#elif defined(CERNLIB_QMCRY)
      CALL UNPACK (LQ(LIN),32,LQ(L4STOX),MAXREX)
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      CALL UPAK32 (LQ(LIN),LQ(L4STOX),MAXREX)
#elif defined(CERNLIB_QMIBMD)
      CALL UPACK (LQ(LIN),LQ(L4STOX),MAXREX)
#elif defined(CERNLIB_QMIRTD)
      CALL UPACK (LQ(LIN),LQ(L4STOX),MAXREX)
#elif 1
      CALL VXINVB (LQ(LIN),24)
#endif

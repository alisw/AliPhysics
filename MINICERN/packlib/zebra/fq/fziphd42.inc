*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:09  mclareni
* Zebra
*
*
*   unpack 60-bit words to 32-bit
*
* fziphd42.inc
*
#if defined(CERNLIB_QMCDC)
      CALL UPKCH (LQ(LIN),LQ(L4STAI),NW4USE,IPKVEC)
*   unpack 64-bit words to 32-bit
*
* fziphd42.inc
*
#elif defined(CERNLIB_QMCRY)
      CALL UNPACK (LQ(LIN),32,LQ(L4STAI),NW4USE)
*   unpack 64-bit words to 32-bit
*
* fziphd42.inc
*
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      CALL UPAK32 (LQ(LIN),LQ(L4STAI),NW4USE)
*   unpack 64-bit words to 32-bit
*
* fziphd42.inc
*
#elif defined(CERNLIB_QMIBMD)
      CALL UPACK (LQ(LIN),LQ(L4STAI),NW4USE)
*   unpack 64-bit words to 32-bit
*
* fziphd42.inc
*
#elif defined(CERNLIB_QMIRTD)
      CALL UPACK (LQ(LIN),LQ(L4STAI),NW4USE)
#elif 1
      CALL VXINVB (LQ(LIN),NW4USE)
#endif

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
* fziphm42.inc
*
#if defined(CERNLIB_QMCDC)
      IF (IUPAKI.NE.0)             GO TO 44
      CALL UPKCH (LQ(LBUF),LQ(L4STAI),NW4IN,IPKVEC)
      GO TO 46
#elif defined(CERNLIB_QMCRY)
      IF (IUPAKI.NE.0)             GO TO 44
      CALL UNPACK (LQ(LBUF),32,LQ(L4STAI),NW4IN)
      GO TO 46
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      IF (IUPAKI.NE.0)             GO TO 44
      CALL UPAK32 (LQ(LBUF),LQ(L4STAI),NW4IN)
      GO TO 46
#elif defined(CERNLIB_QMIBMD)
      IF (IUPAKI.NE.0)             GO TO 44
      CALL UPACK (LQ(LBUF),LQ(L4STAI),NW4IN)
      GO TO 46
#elif defined(CERNLIB_QMIRTD)
      IF (IUPAKI.NE.0)             GO TO 44
      CALL UPACK (LQ(LBUF),LQ(L4STAI),NW4IN)
      GO TO 46
#elif 1
      IF (IUPAKI.NE.0)             GO TO 44
      CALL VXINVC (LQ(LBUF),LQ(L4STAI),NW4IN)
      GO TO 46
#endif

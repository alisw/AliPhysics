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
* fzophr22.inc
*
#if defined(CERNLIB_QMCDC)
      IF (NWFILX.NE.0)  NWU = NWR - (32*NWFILX)/60
      CALL PKCHAR (LQ(LWR),LQ(LWR-1),NWPK,IPKVEC)
      LWR  = LWR - 1
*   pack 32-bit words to 64-bits
*
* fzophr22.inc
*
#elif defined(CERNLIB_QMCRY)
      IF (NWFILX.NE.0)  THEN
          NWU  = (NWPK+1) / 2
          NWPK =  NWU * 2
        ENDIF
      CALL PACK (LQ(LWR),32,LQ(LWR),NWPK)
*   pack 32-bit words to 64-bits
*
* fzophr22.inc
*
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
      IF (NWFILX.NE.0)  THEN
          NWU  = (NWPK+1) / 2
          NWPK =  NWU * 2
        ENDIF
      CALL PACK32 (LQ(LWR),LQ(LWR),NWPK)
*   pack 32-bit words to 64-bits
*
* fzophr22.inc
*
#elif defined(CERNLIB_QMIBMD)
      IF (NWFILX.NE.0)  THEN
          NWU  = (NWPK+1) / 2
          NWPK =  NWU * 2
        ENDIF
      CALL PACK (LQ(LWR),LQ(LWR),NWPK)
*   pack 32-bit words to 64-bits
*
* fzophr22.inc
*
#elif defined(CERNLIB_QMIRTD)
      IF (NWFILX.NE.0)  THEN
          NWU  = (NWPK+1) / 2
          NWPK =  NWU * 2
        ENDIF
      CALL PACK (LQ(LWR),LQ(LWR),NWPK)
#elif 1
      CALL VXINVB (LQ(LWR),NWPK)
#endif

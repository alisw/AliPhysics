*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*
* fzfilelg.inc
*
#if defined(CERNLIB_QMIRTD)
      IF (IUPAK.EQ.0)  THEN
          NWRECX = (32*NWRECX) / 64
          JOFFSI = LRECX+2 - NWRECX
        ENDIF
#else
      IF (IUPAK.EQ.0)  THEN
          NWRECX = (32*NWRECX) / IQBITW
          JOFFSI = LRECX+2 - NWRECX
        ENDIF
#endif

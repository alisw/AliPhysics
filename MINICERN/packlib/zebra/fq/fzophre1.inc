*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzophre1.inc
*
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
#if !defined(CERNLIB_QIBM_VM)
      IF (IACMOX.EQ.0)             GO TO 87
      CALL IODD (LUNX,ISW)
      IQUEST(11) = 1
      GO TO 89
#endif
#endif
#if (defined(CERNLIB_QMIBM))&&(defined(CERNLIB_IOPACK))
      IF (IACMOX.EQ.0)             GO TO 87
      CALL IOMOVE (LUNX,-1,ISW)
      IQUEST(11) = 1
      GO TO 89
*    --------------  sequences for input   ---------------------
#endif
*
* fzophre1.inc
*
#if defined(CERNLIB_QMND3)
      IF (IACMOX.EQ.0)             GO TO 87
      DO 84 J=1,NEOFU
      NIL = 1
      CALL MAGTP (12B,NIL,LUNX,NIL,NIL)
   84 CONTINUE
      GO TO 89
*    --------------  conversions for output  -------------------
#endif

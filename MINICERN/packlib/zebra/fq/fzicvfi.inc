*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:09  mclareni
* Zebra
*
*
*   sign-extension for integer input
*
* fzicvfi.inc
*
#if defined(CERNLIB_QMCDC)
  201 DO 229 J=1,NWDO
      ITHA = MS(JMS+1)
      IF ((ITHA.AND.O"20000000000").NE.0)  ITHA= MASK(28).OR.(ITHA-1)
      MT(JMT+1) = ITHA
      JMT = JMT + 1
  229 JMS = JMS + 1
      GO TO 801
#endif
*    use default CALL IE3TOS for input single prec.
#ifndef CERNLIB_FQIE3TSC
#endif
*    use default CALL IE3TOD for input double prec.
#ifndef CERNLIB_FQIE3TDC
#endif
*   sign-extension for integer input
*
* fzicvfi.inc
*
#if defined(CERNLIB_QMCRY)
  201 DO 229 J=1,NWDO
      ITHA = MS(JMS+1)
      IF ((ITHA.AND.20000000000B).NE.0)  ITHA= ITHA .OR. MASK(32)
      MT(JMT+1) = ITHA
      JMT = JMT + 1
  229 JMS = JMS + 1
      GO TO 801
#endif
*     use default in-line conversion for input floating
#ifndef CERNLIB_FQCVDFF
#endif
*   sign-extension for integer input
*
* fzicvfi.inc
*
#if (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
  201 DO 229 J=1,NWDO
      ITHA = MS(JMS+1)
      IF ((ITHA.AND.MSKB32).NE.0)  ITHA= ITHA .OR. MSKL32
      MT(JMT+1) = ITHA
      JMT = JMT + 1
  229 JMS = JMS + 1
      GO TO 801
#endif
*    use default CALL IE3TOS for input single prec.
#ifndef CERNLIB_FQIE3TSC
#endif
*    use default CALL IE3TOD for input double prec.
#ifndef CERNLIB_FQIE3TDC
#endif

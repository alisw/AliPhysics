*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:08  mclareni
* Zebra
*
*
*
* fzocvfi.inc
*
#if defined(CERNLIB_QMCDC)
C--      CDC :  two's complement for integer output
  201 JPRO = 0
      DO 229 J=1,NWDO
      ITHA = MS(JMS+1)
      ITHB = MASK(28) .AND. ITHA
      IF (ITHB.EQ.0)                  GO TO 227
      IF ((ITHB.XOR.MASK(28)).EQ.0)   GO TO 226
      JPRO = JMS
      GO TO 227

  226 ITHA = ((.NOT.MASK(28)).AND.ITHA) + 1
  227 MT(JMT+1) = ITHA
      JMT = JMT + 1
  229 JMS = JMS + 1
      IF (JPRO.EQ.0)               GO TO 801
      IFOCON(1) = 2
      IFOCON(2) = JPRO
      IFOCON(3) = MS(JPRO+1)
      GO TO 801
#endif
*    use default CALL IE3FOS for output single prec.
#ifndef CERNLIB_FQIE3FSC
#endif
*    use default CALL IE3FOD for output double prec.
#ifndef CERNLIB_FQIE3FDC
#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:12  mclareni
* Zebra
*
*
*   cv IEEE <- single
*
* fzocvff.inc
*
#if defined(CERNLIB_FQIE3FSC)
C--       Default conversion from internal to single IEEE
  301 CALL IE3FOS (MS(JMS+1),MT(JMT+1),NWDO,JBAD)
      IF (JBAD.NE.0)               GO TO 331
  304 JMT = JMT + NWDO
      JMS = JMS + NWDO
      GO TO 801

  331 JBAD = JBAD-1 + JMS
      IFOCON(1) = 3
      IFOCON(2) = JBAD
      IFOCON(3) = MS(JBAD+1)
      GO TO 304
#endif

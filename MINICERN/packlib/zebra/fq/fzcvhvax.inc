*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*   invert ASCII characters
*
* fzcvhvax.inc
*
#if defined(CERNLIB_QMVAX)
C--       VAX : invert ASCII characters
      CALL VXINVC (MS(JMS+1),MT(JMT+1),NWDO)
      JMT = JMT + NWDO
      JMS = JMS + NWDO
      GO TO 801
#endif

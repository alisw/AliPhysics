*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*   swop halves of double precision numbers
*
* fzcvdtmo.inc
*
#if defined(CERNLIB_QMTMO)
C--       TMO : exchange left and right half of each d/p number
      DO 449  JL=1,NDPN
      MT(JMT+1) = MS(JMS+2)
      MT(JMT+2) = MS(JMS+1)
      JMT = JMT + 2
  449 JMS = JMS + 2
#endif

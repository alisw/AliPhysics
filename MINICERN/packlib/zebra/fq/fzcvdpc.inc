*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:19:49  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:39  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   swop halves of double precision numbers
*
* fzcvdpc.inc
*
#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
C--       PC : exchange left and right half of each d/p number
      DO 449  JL=1,NDPN
      MT(JMT+1) = MS(JMS+2)
      MT(JMT+2) = MS(JMS+1)
      JMT = JMT + 2
  449 JMS = JMS + 2
#endif

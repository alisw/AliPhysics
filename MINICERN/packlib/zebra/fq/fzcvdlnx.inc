*
* $Id$
*
* $Log$
* Revision 1.2  1998/09/25 09:32:54  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   swop halves of double precision numbers
*
* fzcvdlnx.inc
*
#if (defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))
C--       LNX : exchange left and right half of each d/p number
      DO 449  JL=1,NDPN
      MT(JMT+1) = MS(JMS+2)
      MT(JMT+2) = MS(JMS+1)
      JMT = JMT + 2
  449 JMS = JMS + 2
#endif

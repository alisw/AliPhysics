*
* $Id$
*
* $Log$
* Revision 1.2  1998/09/25 09:32:56  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.1.1.1  1996/03/06 10:47:11  mclareni
* Zebra
*
*
*   invert ASCII characters
*
* fzcvhlnx.inc
*
#if (defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))
C--       VAX : invert ASCII characters
      CALL VXINVC (MS(JMS+1),MT(JMT+1),NWDO)
      JMT = JMT + NWDO
      JMS = JMS + NWDO
      GO TO 801
#endif

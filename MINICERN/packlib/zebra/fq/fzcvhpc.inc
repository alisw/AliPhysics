*
* $Id$
*
* $Log$
* Revision 1.2  1997/03/14 17:19:50  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:40  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   invert ASCII characters
*
* fzcvhpc.inc
*
#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
C--       PC : invert ASCII characters
      CALL VXINVC (MS(JMS+1),MT(JMT+1),NWDO)
      JMT = JMT + NWDO
      JMS = JMS + NWDO
      GO TO 801
#endif

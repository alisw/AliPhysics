*
* $Id$
*
* $Log$
* Revision 1.3  1998/09/25 09:33:04  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.2  1997/03/14 17:19:51  mclareni
* WNT mods
*
* Revision 1.1.1.1.2.1  1997/01/21 11:33:42  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*                cv IEEE <- double
*
* fzocvfd.inc
*
#if defined(CERNLIB_QMDOS)||(defined(CERNLIB_QMLNX) && !defined(CERNLIB_PPC))||defined(CERNLIB_QMTMO)||defined(CERNLIB_QMVMI) || defined(CERNLIB_WINNT)
C--       : exchange left and right half of each d/p number
      DO 449  JL=1,NDPN
      MT(JMT+1) = MS(JMS+2)
      MT(JMT+2) = MS(JMS+1)
      JMT = JMT + 2
  449 JMS = JMS + 2
#elif defined(CERNLIB_FQIE3FDC)
C--       Default conversion from internal to double IEEE
      CALL IE3FOD (MS(JMS+1),MT(JMT+1),NDPN,JBAD)
      IF (JBAD.NE.0)  THEN
          JBAD = 2*JBAD-2 + JMS
          IFOCON(1) = 4
          IFOCON(2) = JBAD
          IFOCON(3) = MS(JBAD+1)
        ENDIF
      JMT = JMT + NWDODB
      JMS = JMS + NWDODB
#elif 1
      CALL UCOPY (MS(JMS+1),MT(JMT+1),NWDODB)
      JMS = JMS + NWDODB
      JMT = JMT + NWDODB
#endif

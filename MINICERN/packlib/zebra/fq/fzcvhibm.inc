*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*   convert EBCDIC <-> ASCII
*
* fzcvhibm.inc
*
#if defined(CERNLIB_QMIBM)
C--       IBM conversion  EBCDIC <--> ASCII
      NTHR = NWDO
  522 NTH  = MIN (NTHR,24)
      NCH  = 4*NTH
      CALL UPKCH (MS(JMS+1),IQUEST,NCH,IPKVEC)
      DO 529  JCH=1,NCH
      J = MIN (MAX(IQUEST(JCH)-JCHMIN,0), JCHMAX)
  529 IQUEST(JCH) = MTAB(J+1)
      CALL PKCHAR (IQUEST,MT(JMT+1),NCH,IPKVEC)
      JMT  = JMT + NTH
      JMS  = JMS + NTH
      NTHR = NTHR - NTH
      IF (NTHR.NE.0)               GO TO 522
      GO TO 801
#endif
*   convert EBCDIC <-> ASCII
*
* fzcvhibm.inc
*
#if defined(CERNLIB_QMIBMD)
C--       IBM conversion  EBCDIC <--> ASCII
      NTHR = NWDO
  522 NTH  = MIN (NTHR,24)
      NCH  = 4*NTH
      CALL UPKCH (MS(JMS+1),IQUEST,NCH,IPKVEC)
      DO 529  JCH=1,NCH
      J = MIN (MAX(IQUEST(JCH)-JCHMIN,0), JCHMAX)
  529 IQUEST(JCH) = MTAB(J+1)
      CALL PKCHAR (IQUEST,MT(JMT+1),NCH,IPKVEC)
      JMT  = JMT + NTH
      JMS  = JMS + NTH
      NTHR = NTHR - NTH
      IF (NTHR.NE.0)               GO TO 522
      GO TO 801
#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzocvfai.inc
*
#if defined(CERNLIB_QMCRY)
C--       Cray B + F : copy as is, clear left half-word
C-        First, check bits lost
C     DO 104  J=1,NWDO
C     IF (SHIFTR(MS(JMS+J),32).NE.0)  GO TO 105
C 104 CONTINUE
C     GO TO 106
C
C-        Conversion problem
C 105 IFOCON(1) = ITYPE
C     IFOCON(2) = JMS + J
C     IFOCON(3) = MS(JMS+J)

C-        Copy right halves only
  106 DO 107  J=1,NWDO
  107 MT(JMT+J) = MS(JMS+J) .AND. 37777777777B
      JMT = JMT + NWDO
      JMS = JMS + NWDO
#elif (defined(CERNLIB_QMCV64))&&(!defined(CERNLIB_QIEEE))
C-        Copy right halves only
  106 DO 107  J=1,NWDO
  107 MT(JMT+J) = MS(JMS+J) .AND. 'FFFFFFFF'X
      JMT = JMT + NWDO
      JMS = JMS + NWDO
#elif 1
      CALL UCOPY (MS(JMS+1),MT(JMT+1),NWDO)
      JMT = JMT + NWDO
      JMS = JMS + NWDO
#endif

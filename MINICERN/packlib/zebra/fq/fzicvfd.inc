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
*                cv IEEE -> double
*
* fzicvfd.inc
*
#if defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
#include "fzcvdpc.inc"
#elif defined(CERNLIB_QMLNX)
#include "fzcvdlnx.inc"
#elif defined(CERNLIB_QMTMO)
#include "fzcvdtmo.inc"
#elif defined(CERNLIB_QMVMI)
#include "fzcvdvmi.inc"
#elif defined(CERNLIB_FQIE3TDC)
C--       Default conversion from double IEEE to internal
      CALL IE3TOD (MS(JMS+1),MT(JMT+1),NDPN,JBAD)
      IF (JBAD.NE.0)  THEN
          JBAD = 2*JBAD - 1
          IFOCON(1) = 4
          IFOCON(2) = JMT + JBAD
          IFOCON(3) = MS(JMS+JBAD)
        ENDIF
      JMT = JMT + NWDODB
      JMS = JMS + NWDODB
#elif defined(CERNLIB_FQCVDFF)
C--       Default conversion from double IEEE to internal
      DO 449  JL=1,NWDODB,2
      ITHA = MS(JMS+1)
      IF (ITHA.EQ.0)               GO TO 442
      ITHB  = MS(JMS+2)
      JSIGN = JBIT (ITHA,32)
      JEXP  = JBYT (ITHA,21,11)
      JMANT = JBYT (ITHA, 1,20)
      IF (JEXP.EQ.0)               GO TO 431
      IF (JEXP.EQ.2047)            GO TO 433

      JMANT = JMANT + 1 048 576
      JEXP  = JEXP  - 1023
      IF (JEXP.GT.JEXMAX)          GO TO 432
  424 IF (JEXP.LT.JEXMIN)          GO TO 441

      JMANS = JBYT (ITHB,17,16)
      JMANR = JBYT (ITHB, 1,16)

      JPRE  = ISIGN (48,JEXP)
      THDB  = (  DBLE(JMANT)
     +         + DBLE(FLOAT(JMANS)*2.**(-16))
     +         + DBLE(FLOAT(JMANR)*2.**(-32)) ) * 2.**(JPRE-20)
      THDB  = THDB * 2.**(JEXP-JPRE)
      IF (JSIGN.EQ.0)              GO TO 447
      THDB = -THDB
      GO TO 447

  431 IF (JMANT.EQ.0)              GO TO 441
      JEXP = -1022
      GO TO 424

  432 JMANT = 0
  433 IF (JMANT.EQ.0)    THEN
          IF (JSIGN.EQ.0)  THEN
              ITHB = IOVPMD
            ELSE
              ITHB = IOVNMD
            ENDIF
        ELSE
          ITHB = IBADMD
        ENDIF
      IFOCON(1) = 4
      IFOCON(2) = JMT
      IFOCON(3) = ITHA
      ITHA = ITHB
      ITHB = IBADME
      GO TO 447

  441 ITHA = 0
  442 ITHB = 0
  447 MT(JMT+1) = ITHA
      MT(JMT+2) = ITHB
      JMT = JMT + 2
  449 JMS = JMS + 2
#elif 1
      CALL UCOPY (MS(JMS+1),MT(JMT+1),NWDODB)
      JMT = JMT + NWDODB
      JMS = JMS + NWDODB
#endif

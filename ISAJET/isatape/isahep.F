#include "isajet/pilot.h"
      SUBROUTINE ISAHEP(MCONV)
C
C...Purpose: to convert ISAJET event record contents to or from
C...the standard event record common block.
C
C   Thanks to Lynn Garren, Fermilab.
C
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/hepevt.inc"
C
C...for event number
#include "isajet/idrun.inc"
C...initial partons
#include "isajet/pjets.inc"
#include "isajet/primar.inc"
C...partons created during decay
#include "isajet/jetset.inc"
C...particles created in the decay, including final state particles
#include "isajet/partcl.inc"
C
      INTEGER MCONV
      INTEGER ITRANS
      INTEGER I1,IHP,MPART,JET,NPOFF,NJHEP,NWHEP,IMO,IJT
      INTEGER JPMO(2,MXJSET),JPDA(2,MXJSET),JMX(MXJSET),JMN(MXJSET)
      INTEGER JTMO(2,MXPTCL),JTDA(2,MXPTCL)
      INTEGER IP,IJ3,IJ2,IJ1,NSUM2,NSUM1,IPT,I,J,KST,KND,K
C
C
C...Conversion from ISAJET to standard
C
      IF(MCONV.EQ.1) THEN
        NEVHEP = IEVT
C...initial jets
        NHEP = 0
C... W or Z
        IF(IDENTW.NE.0)THEN
          NHEP = NHEP + 1
          ISTHEP(NHEP)=12
          JMOHEP(1,NHEP)=0
          JMOHEP(2,NHEP)=0
          JDAHEP(1,NHEP)= 2
          JDAHEP(2,NHEP)= NJET + 1
          IDHEP(NHEP) = ITRANS(IDENTW,1)
          DO 100 J=1,5
 100        PHEP(J,NHEP) = QWJET(J)
        ENDIF
        NWHEP = NHEP
C... jets
        IF(NJET.GT.0)THEN
          DO 120 I=1,NJET
            NHEP = NHEP + 1
            ISTHEP(NHEP)=11
            JMOHEP(1,NHEP)=0
            IF(IDENTW.NE.0) JMOHEP(1,NHEP) = 1
            JMOHEP(2,NHEP)= I
            JDAHEP(1,NHEP)=0
            JDAHEP(2,NHEP)=0
            IDHEP(NHEP) = ITRANS(IDJETS(I),1)
            DO 110 J=1,5
 110          PHEP(J,NHEP) = PJETS(J,I)
 120      CONTINUE
        ENDIF
        NJHEP = NHEP
C... pairs
        IF(NPAIR.GT.0)THEN
          DO 150 I=1,NPAIR
            NHEP = NHEP + 1
            ISTHEP(NHEP)=13
            JMOHEP(1,NHEP)= JPAIR(I) + NWHEP
            JMOHEP(2,NHEP)= JPAIR(I)
            JDAHEP(1,NHEP)=0
            JDAHEP(2,NHEP)=0
            IDHEP(NHEP) = ITRANS(IDPAIR(I),1)
            DO 140 J=1,5
 140          PHEP(J,NHEP) = PPAIR(J,I)
 150      CONTINUE
        ENDIF
        DO 160 I=1,NHEP
          DO 160 J=1,4
 160        VHEP(J,I) = 0.
C...save offset into hep list
        NPOFF = NHEP
C...partons
        DO 200 I=1,NJSET
          IHP = NHEP + I
C...use JMX and JMN to find daughters in hadron list
          JMX(I) = 0
          JMN(I) = NHEP + NPTCL + 1
          IDHEP(IHP) = ITRANS(JTYPE(I),1)
          MPART=MOD(JORIG(I),JPACK)
          JMOHEP(1,IHP)=0
          IJT = JORIG(I)/JPACK
          IF(MPART.NE.0)THEN
            JMOHEP(1,IHP)=MPART+NHEP
          ELSEIF(MPART.EQ.0 .AND. IJT.LT.10)THEN
C...find mother in jet/pair list
            IMO = IJT + NWHEP
            IF(NJHEP.LT.NPOFF)THEN
              KST = NJHEP + 1
              DO 170 K=KST,NPOFF
                IF(IDHEP(K).EQ.IDHEP(IHP)) IMO=K
 170          CONTINUE
            ENDIF
            JMOHEP(1,IHP)= IMO
            IF(JDAHEP(1,IMO).EQ.0) JDAHEP(1,IMO)=IHP
            JDAHEP(1,IMO) = MIN(IHP,JDAHEP(1,IMO))
            JDAHEP(2,IMO) = MAX(IHP,JDAHEP(2,IMO))
C...amend information if a parton thinks this is it's daughter
            KND = IHP-1
            DO 175 K=NPOFF,KND
              IF(IHP.GE.JDAHEP(1,K) .AND. IHP.LE.JDAHEP(2,K))
     1             JMOHEP(1,IHP)=K
 175        CONTINUE
          ENDIF
          JMOHEP(2,IHP)= IJT
          IF(JDCAY(I).EQ.0)THEN
            ISTHEP(IHP) = 21
            JDAHEP(1,IHP)=0
            JDAHEP(2,IHP)=0
          ELSE
            ISTHEP(IHP) = 22
            JDAHEP(1,IHP) = JDCAY(I)/JPACK + NHEP
            JDAHEP(2,IHP) = MOD(JDCAY(I),JPACK) + NHEP
          ENDIF
          DO 180 J=1,5
 180        PHEP(J,IHP) = PJSET(J,I)
          DO 190 J=1,4
 190        VHEP(J,IHP) = 0.
 200    CONTINUE
        NHEP = NHEP + NJSET
C...hadrons
        DO 250 I=1,NPTCL
          IHP = NHEP + I
          IDHEP(IHP) = ITRANS(IDENT(I),1)
          I1 = MOD(IABS(IORIG(I)),IPACK)
          JMOHEP(1,IHP)=0
          JMOHEP(2,IHP)=IABS(IORIG(I))/IPACK
C...mother is pomeron
          IF(I1.EQ.0)THEN
C...mother is in parton list
          ELSEIF(IORIG(I).LT.0)THEN
            JMOHEP(1,IHP) = I1 + NPOFF
            JMN(I1) = MIN(JMN(I1),I)
            JMX(I1) = MAX(JMX(I1),I)
C...mother is in hadron list
          ELSEIF(IORIG(I).GT.0)THEN
            JMOHEP(1,IHP) = I1 + NHEP
          ENDIF
          IF(IDCAY(I).EQ.0)THEN
            ISTHEP(IHP) = 1
            JDAHEP(1,IHP)=0
            JDAHEP(2,IHP)=0
          ELSE
            ISTHEP(IHP) = 2
            JDAHEP(1,IHP) = IDCAY(I)/IPACK + NHEP
            JDAHEP(2,IHP) = MOD(IDCAY(I),IPACK) + NHEP
          ENDIF
          DO 210 J=1,5
 210        PHEP(J,IHP) = PPTCL(J,I)
          DO 220 J=1,4
 220        VHEP(J,IHP) = 0.
 250    CONTINUE
        NHEP = NHEP + NPTCL
C...fill in missing daughter info for partons
        DO 270 I=1,NJSET
          IF(JMX(I).NE.0)THEN
            JDAHEP(1,I+NPOFF) = JMN(I) + NPOFF + NJSET
            JDAHEP(2,I+NPOFF) = JMX(I) + NPOFF + NJSET
          ENDIF
 270    CONTINUE
C
C...Conversion from standard to ISAJET
C
      ELSEIF(MCONV.EQ.2)THEN
        IEVT = NEVHEP
C...        missing information
        IDENTW = 0
        NPAIR = 0
        DO 330 I=1,5
          QWJET(I) = 0.
          DO 330 J=1,4
            PPAIR(I,J) = 0.
 330    CONTINUE
        DO 340 I=1,4
          IDPAIR(I) = 0
 340      JPAIR(I) = 0
C...zero counters
        IJ1 = 0
        IJ2 = 0
        IJ3 = 0
        IP = 0
        IPT = 0
        DO 500 I=1,NHEP
C...initial jets
C... jets
          IF(ISTHEP(I).EQ.11)THEN
            IJ1 = IJ1 + 1
            IDJETS(IJ1) = ITRANS(IDHEP(I),2)
            DO 410 J=1,5
 410          PJETS(J,IJ1) = PHEP(J,I)
C... W
          ELSEIF(ISTHEP(I).EQ.12)THEN
            IJ2 = IJ2 + 1
            IDENTW = ITRANS(IDHEP(I),2)
            DO 420 J=1,5
 420          QWJET(J) = PHEP(J,I)
C... pairs
          ELSEIF(ISTHEP(I).EQ.13)THEN
            IJ3 = IJ3 + 1
            IDPAIR(IJ3) = ITRANS(IDHEP(I),2)
            JPAIR(IJ3) = JMOHEP(2,I)
            DO 430 J=1,5
 430          PPAIR(J,IJ3) = PHEP(J,I)
C...partons
          ELSEIF(ISTHEP(I).EQ.21 .OR. ISTHEP(I).EQ.22)THEN
            IP = IP + 1
            JTYPE(IP) = ITRANS(IDHEP(I),2)
            DO 440 J=1,5
 440          PJSET(J,IP) = PHEP(J,I)
C...      temporary storage until have counts
            JPMO(1,IP) = JMOHEP(1,I)
            JPMO(2,IP) = JMOHEP(2,I)
            JPDA(1,IP) = JDAHEP(1,I)
            JPDA(2,IP) = JDAHEP(2,I)
C...hadrons
          ELSE
            IPT = IPT + 1
            IDENT(IPT) = ITRANS(IDHEP(I),2)
            DO 450 J=1,5
 450          PPTCL(J,IPT) = PHEP(J,I)
C...      temporary storage until have counts
            JTMO(1,IPT) = JMOHEP(1,I)
            JTMO(2,IPT) = JMOHEP(2,I)
            JTDA(1,IPT) = JDAHEP(1,I)
            JTDA(2,IPT) = JDAHEP(2,I)
          ENDIF
 500    CONTINUE
C...completed counts
        NJET = IJ1
        NPAIR = IJ3
        NJSET = IP
        NPTCL = IPT
C...get mother/daughter information
        NSUM1 = NJET + IJ2 + NPAIR
        NSUM2 = NSUM1 + NJSET
        DO 520 I=1,NJSET
          IF(JPDA(1,I).EQ.0)THEN
            JDCAY(I) = 0
          ELSEIF(JPDA(1,I).GT.NSUM2)THEN
            JDCAY(I) = 0
          ELSE
            JDCAY(I) = JPACK*(JPDA(1,I)-NSUM1) + JPDA(2,I)-NSUM1
          ENDIF
          IF(JPMO(1,I).LE.NSUM1)THEN
            JORIG(I) = JPACK*JPMO(2,I)
          ELSE
            JORIG(I) = JPACK*JPMO(2,I) + JPMO(1,I)-NSUM1
          ENDIF
 520    CONTINUE
        DO 550 I=1,NPTCL
          IF(JTDA(1,I).EQ.0)THEN
            IDCAY(I) = 0
          ELSE
            IDCAY(I) = IPACK*(JTDA(1,I)-NSUM2) + JTDA(2,I)-NSUM2
          ENDIF
          IF(JTMO(1,I).LE.NSUM1)THEN
            IORIG(I) = JTMO(2,I)*IPACK + 0
          ELSEIF(JTMO(1,I).LE.NSUM2)THEN
            IORIG(I) = -(JTMO(2,I)*IPACK + JTMO(1,I)-NSUM1)
          ELSE
            IORIG(I) = JTMO(2,I)*IPACK + JTMO(1,I)-NSUM2
          ENDIF
 550    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE PHOBOS(IP,PBOOS1,PBOOS2,FIRST,LAST)
C.----------------------------------------------------------------------
C.
C.    PHOBOS:   PHOton radiation in decays BOoSt routine
C.
C.    Purpose:  Boost particles  in  cascade decay  to parent rest frame
C.              and boost back with modified boost vector.
C.
C.    Input Parameters:       IP:  pointer of particle starting chain
C.                                 to be boosted
C.                        PBOOS1:  Boost vector to rest frame,
C.                        PBOOS2:  Boost vector to modified frame,
C.                        FIRST:   Pointer to first particle to be boos-
C.                                 ted (/PH_HEPEVT/),
C.                        LAST:    Pointer to last  particle to be boos-
C.                                 ted (/PH_HEPEVT/).
C.
C.    Output Parameters:  Common /PH_HEPEVT/.
C.
C.    Author(s):  B. van Eijk                     Created at:  13/02/90
C.                Z. Was                          Last Update: 16/11/93
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION BET1(3),BET2(3),GAM1,GAM2,PB,DATA
      INTEGER I,J,FIRST,LAST,MAXSTA,NSTACK,IP
      PARAMETER (MAXSTA=10000)
      INTEGER STACK(MAXSTA)
      REAL*8 PBOOS1(5),PBOOS2(5)
      INTEGER NMXHEP
      PARAMETER (NMXHEP=10000)
      INTEGER IDHEP,ISTHEP,JDAHEP,JMOHEP,NEVHEP,NHEP
      REAL*8 PHEP,VHEP
      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      IF ((LAST.EQ.0).OR.(LAST.LT.FIRST)) RETURN
      NSTACK=0
      DO 10 J=1,3
        BET1(J)=-PBOOS1(J)/PBOOS1(5)
   10 BET2(J)=PBOOS2(J)/PBOOS2(5)
      GAM1=PBOOS1(4)/PBOOS1(5)
      GAM2=PBOOS2(4)/PBOOS2(5)
C--
C--   Boost vector to parent rest frame...
   20 DO 50 I=FIRST,LAST
        PB=BET1(1)*PHEP(1,I)+BET1(2)*PHEP(2,I)+BET1(3)*PHEP(3,I)
        IF (JMOHEP(1,I).EQ.IP) THEN
         DO 30 J=1,3
   30    PHEP(J,I)=PHEP(J,I)+BET1(J)*(PHEP(4,I)+PB/(GAM1+1.D0))
         PHEP(4,I)=GAM1*PHEP(4,I)+PB
C--
C--    ...and boost back to modified parent frame.
         PB=BET2(1)*PHEP(1,I)+BET2(2)*PHEP(2,I)+BET2(3)*PHEP(3,I)
         DO 40 J=1,3
   40    PHEP(J,I)=PHEP(J,I)+BET2(J)*(PHEP(4,I)+PB/(GAM2+1.D0))
         PHEP(4,I)=GAM2*PHEP(4,I)+PB
         IF (JDAHEP(1,I).NE.0) THEN
           NSTACK=NSTACK+1
C--
C--    Check on stack length...
           IF (NSTACK.GT.MAXSTA) THEN
             DATA=NSTACK
             CALL PHOERR(7,'PHOBOS',DATA)
           ENDIF
           STACK(NSTACK)=I
         ENDIF
        ENDIF
   50 CONTINUE
      IF (NSTACK.NE.0) THEN
C--
C--   Now go one step further in the decay tree...
        FIRST=JDAHEP(1,STACK(NSTACK))
        LAST=JDAHEP(2,STACK(NSTACK))
        IP=STACK(NSTACK)
        NSTACK=NSTACK-1
        GOTO 20
      ENDIF
      RETURN
      END

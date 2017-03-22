*CMZ :          17/07/95  11.18.11  by  Jonathan Butterworth
*-- Author : JMB  03/12/92
      SUBROUTINE HWNSCT(J)
C --------------------------------------------------------------------------
C Purpose: Choose  how many scatters in this event (NSCAT).
C          If the number of scatters returned is 0, the event
C          should be rejected due to the reductiion of the cross
C          by eikonalisation.
C          J (return argument) is the index of the upper edge of the
C          bin in s in which we lie.
C          If the chosen particles are at a fixed cms, J is always 1
C --------------------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"

      DOUBLE PRECISION Z, EIK, UNEIK, JMGAMZ, TEST
      DOUBLE PRECISION JMRNDM, DICE, SOFAR, PROBS(MAXMS)
      INTEGER I, J, DUM, IERR
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST
      REAL tmp1, tmp2

      NSCAT=0
      IF (IERROR.NE.0) RETURN

      IERR = 0

C --  Z is the fraction of the beam-beam centre-of-mass energy (squared)
C     available in the subcollision. 
      IF (JCMVAR.EQ.0) THEN
        Z=-1.D0
        J=1
      ELSE
        Z=PHEP(5,3)**2/(2.D0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2))
      ENDIF
C
C --  Stage one: Keep the event?
      IF (JMUEO.NE.0) THEN
         TEST = 1.
      ELSE
         EIK=JMGAMZ(Z,2,J)
         UNEIK=JMGAMZ(Z,1,J)
         TEST=EIK/UNEIK      
C     Adjust the probability of zero scatters to the fact that some
C     scatters are lost in event generation.
         IF (JMARRY(5+MAXMS,J).GT.0.AND.JMUEO.EQ.0) THEN
            TEST=TEST*JMARRY(6+MAXMS,J)/JMARRY(5+MAXMS,J)
         ENDIF
         DICE = JMRNDM(2,DUM)
         
         IF (DICE.GT.TEST) THEN
C     Event is rejected
            NSCAT = 0
            RETURN
         ENDIF
      ENDIF

C --  Stage two: Decide how many scatters in this event.
C     First find the respective probabilities.
      CALL JMPN( PROBS, Z, J )
      DICE = JMRNDM(2,DUM)
      SOFAR = 0.D0
      IF (JMBUG.GT.0) THEN
         IF (FIRST) THEN
            tmp1=0.
            tmp2=0.
            do i=1,maxms
               write(*,*) 'HWNSCT:Prob ',i,'=',probs(i)
               tmp1=tmp1+probs(i)
               tmp2=tmp2+probs(i)*float(i)
            enddo
            write(*,*) 'Summed prob =',tmp1
            write(*,*) 'Weighted average prob =',tmp2/float(maxms)
            first=.false.
         ENDIF
         IF (JMBUG.GT.2) write(*,*) 'HWNSCT: dice=',dice
      ENDIF
      DO I = 1, MAXMS
        SOFAR = SOFAR + PROBS(I)
        IF (DICE.LT.SOFAR) THEN
          NSCAT = I
          GOTO 11
        ENDIF
      ENDDO
      WRITE(JMOUT,*)
     &    '*** HWNSCT WARNING:TOO MANY SCATTERS IN EVENT.'
      NSCAT = MAXMS

 11   CONTINUE

      END




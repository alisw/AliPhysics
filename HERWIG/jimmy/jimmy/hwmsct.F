      SUBROUTINE HWMSCT(ABORT)
C -----------------------------------------------------------------
C Administer multiple scattering
C
C If the returned argument (ABORT) is TRUE then the event has been
C vetoed by the eikonalisation, and should be abandoned.
C -----------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"

      LOGICAL ABORT
      INTEGER N, I, J, REPORT, CHECKSUM, CHECKEV, COUNTER
      INTEGER TMPPR
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST, CHECKSUM, CHECKEV

      IF (FIRST) THEN
         CHECKSUM=0.
         CHECKEV=0.
         FIRST=.FALSE.
      ENDIF

C -- Abort the event if there are any errors.
      IF (IERROR.NE.0) THEN
        ABORT=.TRUE.
        CALL HWUFNE
        RETURN
      ENDIF

      ABORT = .FALSE.

C     Decide how many hard scatters
C     If NSCAT=0, this event was rejected (this can happen because the 
C     cross section at a given centre-of-mass energy changes with 
C     eikonalisation).
C     J returns the index of the upper S bound.
      CALL HWNSCT(J)

      IF (JMBUG.GT.0) THEN
         WRITE(*,*) 'HWMSCT:NUMBER OF SCATTERS REQUESTED =',NSCAT    
         CHECKSUM=CHECKSUM+NSCAT
         IF (NSCAT.GT.0) THEN
            CHECKEV=CHECKEV+1
            WRITE(*,*) 'RUNNING AVERAGE=',FLOAT(CHECKSUM)/FLOAT(CHECKEV)
         ENDIF
      ENDIF

      IF (NSCAT.EQ.0) THEN
C --     Finish event (making sure aborted events aren't printed)
         TMPPR=MAXPR
         MAXPR=0
         CALL HWUFNE
         MAXPR=TMPPR
         ABORT = .TRUE.
         NEVHEP=NEVHEP-1
         RETURN
      ENDIF

      REPORT = 0
      COUNTER = 0
      N = NSCAT

      DO I = 2, N

 10      CONTINUE

         IF (REPORT.EQ.0) THEN
            
            CALL HWHSCT(REPORT,(I.EQ.2),JMUEO,PTJIM)

            IF (REPORT.EQ.5) THEN
               WRITE(*,*) 'FATAL ERROR'
               RETURN
            ENDIF

c     If there's an error which isn't a fatal error, retry.
            IF (REPORT.NE.0) THEN

               IF (JMBUG.GT.0) WRITE(*,*) 'Report=',report
               REPORT=0
               IF (COUNTER.LT.MAXMSTRY) THEN
                  COUNTER=COUNTER+1
                  GOTO 10
               ELSE
c     Avoid infinite loops. If there's really no phase space, reduce
c     the number of scatters.     
                  NSCAT=NSCAT-1
                  IF (JMBUG.GT.0) WRITE(*,*) 'LOST A SCATTER'
                  COUNTER=0
               ENDIF
            ENDIF

            IF (ANOMOFF) THEN
               
               IF (ANOMSC(1,1).NE.0.OR.ANOMSC(1,2).NE.0)THEN
                  WRITE(JMOUT,*) 'Anom. scat',ANOMSC(1,1),ANOMSC(1,2)
                  NSCAT=NSCAT-1
                  REPORT=6
               ENDIF
               
            ENDIF
            
         ENDIF
         
      ENDDO

      IF (JMBUG.GT.0) THEN
         WRITE(*,*) 'HWMSCT:NUMBER OF SCATTERS GENERATED =',NSCAT
      ENDIF

 100  CONTINUE

C     =====================================================================
C     Store the number of events lost from the HERWIG cross section.
C     (which equals the number of "multiple" scatters)
      TOTSCAT = TOTSCAT + NSCAT
      IF (JMUEO.EQ.0) THEN
         NLOST   = NLOST + (NSCAT - 1.D0)
      ENDIF
C     Get the S-hat distribution correct. To do this we must alter
C     the amount of vetoing of events, as this is based initially on
C     the simple eikonal model without taking into account "lost"
C     scatters.
C     We should decrease the chances of an event being rejected by
C     a factor of NSCAT(S)/N(S) where NSCAT is the number of scatters
C     actually generated at this S, and N(S) is the number that would have
C     been generated so far by the "simple" model.
C     Store these numbers in JMARRY(MAXMS+5,I) and JMARRY(MAXMS+6,I)

c     Turn this feature off if the process being generated is not
c     eikonalised itself (i.e. MI are being used only for the
c     underlying event.
      JMARRY(MAXMS+5,J)   = JMARRY(MAXMS+5,J)+FLOAT(NSCAT)
      JMARRY(MAXMS+6,J)   = JMARRY(MAXMS+6,J)+FLOAT(N)

      RETURN
      END





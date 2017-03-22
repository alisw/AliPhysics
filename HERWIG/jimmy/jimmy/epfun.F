      DOUBLE PRECISION FUNCTION EPFUN( Y, I )
C     -----------------------------------------------------------------
C     Purpose: To construct the inclusive & eikonalised hadronic cross
C     section at the given Z. The bremstralung factor for photons in not
C     included.
C     -----------------------------------------------------------------
#include "herwig65.inc"
#include "jimmy.inc"

      INTEGER I, J
      DOUBLE PRECISION Y, Z, JMSHAT, BREMFAC, JMSN
      DOUBLE PRECISION EIKON, QCDRES, ARG1, ARG2, JMFWW

      Z = Y
      JMSHAT=Z*2.D0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2)

c     Uneikonalised cross section
      ARG1=QCDRES( Z )

c     Eikonalise the cross section
      ARG2=EIKON( ARG1 )

C --  Store the p(n) values too.
      DO J=1,MAXMS
        IF (ARG2.GT.1.0D-17) THEN
           IF (JMUEO.NE.0) THEN
              JMARRY( 4+J, I ) = JMSN( ARG1, J )/ARG1
           ELSE
              JMARRY( 4+J, I ) = JMSN( ARG1, J )/ARG2
           ENDIF
        ELSE
          JMARRY( 4+J, I ) = 0.D0
        ENDIF
      ENDDO

      IF (I.LE.NPSIMP) THEN

C       Store (and write out) intermediate results
C       for the xsecn.
        WRITE(JMOUT,801) ARG1*GEV2NB/1000.0,ARG2*GEV2NB/1000.0, JMSHAT
        JMARRY( 1, I ) = Z
C       Store uneikonalised cross section.
        JMARRY( 2, I ) = ARG1
C       Store eikonalised cross section.
        JMARRY( 3, I ) = ARG2

      ELSE

        WRITE(JMOUT,*) 'TOO MANY CALLS TO EPFUN!-FATAL'
        STOP

      ENDIF

      EPFUN = ARG2

      IF (JMBUG.GT.2) THEN
c         WRITE(JMOUT,*) Z,ARG1*BREMFAC,ARG2*BREMFAC
         WRITE(JMOUT,*) Z,ARG1,ARG2
     &,JMARRY(5,I),JMARRY(6,I),JMARRY(7,I),JMARRY(8,I),JMARRY(9,I)
      ENDIF

      RETURN
 801  FORMAT(2X,'HADRONIC C-S: UNEIKONALIZED=',G8.2
     &     ,'ub, EIKONALIZED=',G8.2,'ub, at s =',G8.2,' GeV2')
      END



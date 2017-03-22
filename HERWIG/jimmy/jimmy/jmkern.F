*-- Author : JRF/JMB
      DOUBLE PRECISION FUNCTION JMKERN( NDIM, V )
C -------------------------------------------------------------
C   Purpose: The kernel for the integral over PT,X1,X2
C            The integral goes from pt2tmp -> ptmax2.
C            When using this function in the initialsation,
C            pt2tmp == ptmin**2 and ptmax2 = the kinematic limit.
C            Otherwise they are the pt2 of the current
C            hard scatter and of the last, respectively.
C -------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"


      INTEGER NDIM
      DOUBLE PRECISION PT2, PT2TMP, PTMAX2, V(NDIM)
      DOUBLE PRECISION JMXS1, FAC1, FAC2
      DOUBLE PRECISION FAC3, LX1, LX2, LPT2
      DOUBLE PRECISION X1,X2,JMS

c     Beam-beam centre-of-mass energy.
      JMS = 2.D0*(EBEAM1*EBEAM2+PBEAM1*PBEAM2)

      IF (FN_TYPE.NE.101) THEN
        WRITE(JMOUT,8902) 'JMKERN:THIS SHOULD NEVER HAPPEN!'
        RETURN
      ENDIF

      PT2TMP = PTJIM**2
C     (This value is never used
C     - it is always greater than the kinematic limit)
      PTMAX2 = JMS/4.D0 - RMASS(1)**2

C     Define region of integration
      FAC1 = -DLOG( 4.D0*PT2TMP/(YGAMMA*JMS) )
      IF ( FAC1 .LE. 0.D0 ) THEN
        JMKERN = 0.D0
        GOTO 8900
      ENDIF

      IF (JMBUG.GT.2) THEN
        write(JMOUT,8902) 'v(1),ygamma,jms,pt2tmp,fac1'
        write(JMOUT,*) v(1),ygamma,jms,pt2tmp
     &     ,fac1
      ENDIF

      LX1 = DLOG( 4.D0*PT2TMP/(YGAMMA*JMS)) +V(1)*FAC1
      X1 = DEXP( LX1 )
      IF ( X1 .GE. ONE ) THEN
        JMKERN = 0.D0
        GOTO 8900
      ENDIF

      FAC2 = -DLOG( 4.D0*PT2TMP/(YGAMMA*X1*JMS) )
      IF ( FAC2 .LE. 0.D0 ) THEN
        JMKERN = 0.D0
        GOTO 8900
      ENDIF

      LX2 = DLOG( 4.D0*PT2TMP/(YGAMMA*JMS*X1) ) + V(2)*FAC2
      X2 = DEXP( LX2 )
      IF ( X2 .GE. 1.D0 ) THEN
        JMKERN = 0.D0
        GOTO 8900
      ENDIF

      IF (NDIM.EQ.3) THEN
        PTMAX2 = MIN(JMS*YGAMMA*X1*X2/4.D0,PTMAX2)
        FAC3 = DLOG( PTMAX2 ) - DLOG( PT2TMP )
        IF ( FAC3 .LE. SMALL ) THEN
          JMKERN = 0.D0
          GOTO 8900
        ENDIF

        LPT2 = DLOG( PT2TMP ) + V(3)*FAC3
        PT2  = DEXP( LPT2 )
        IF (JMBUG.GT.2) THEN
          write(JMOUT,8902) 'lx1,x1,lx2,x2,pt2'
          write(JMOUT,*) lx1,x1,lx2,x2,pt2
        ENDIF
        JMKERN = X1*X2*PT2*FAC1*FAC2*FAC3*
     &      PHAD * JMXS1(X1,X2,PT2,0,0)

	IF (JMBUG.GT.2) THEN
          write(JMOUT,8902) 'JMKERN,X1,X2,PT2,FAC1,FAC2,FAC3,PHAD'
          write(JMOUT,*) JMKERN,X1,X2,PT2,FAC1,FAC2,FAC3,PHAD
	ENDIF

        IF ( JMKERN .LE. 0.D-30 ) JMKERN = 0.D0

      ELSE

         WRITE(JMOUT,8902) 'JMKERN ERROR'

      ENDIF

 8900 CONTINUE

      IF (JMBUG.GT.2) THEN
        WRITE(JMOUT,8901) JMKERN,X1,X2,PT2
      ENDIF
 8901 FORMAT
     & (1X,'XS=',F12.6,' X1=',F10.8,' X2=',F10.8,' PT2=',F10.3)
 8902 FORMAT(A)
      RETURN
      END

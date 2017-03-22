
*     Include file for JIMMY

      DOUBLE PRECISION YGAMMA, JMZMIN, JMRAD, PTJIM
      DOUBLE PRECISION PHAD, JMU2, JMV2, SMALL, JMARRY
c     JMARRY is the array storing gamma-p xsec at various z, & 
c	max weight for each z
      DOUBLE PRECISION TOTSCAT, NLOST
      INTEGER MAXMS, NPSIMP, MSFLAG, JMPTYP, JCMVAR, NPROC
      LOGICAL ANOMOFF
      PARAMETER( NPROC = 117 )
      PARAMETER( MAXMS  = 100  )  ! Maximum multiple scatters
      PARAMETER( NPSIMP = 16 )  ! No. of Simpson rule (YBJ)
C                                 intervals (must be even)
      PARAMETER( SMALL  = 1.0D-20  )
      INTEGER JMOUT, JMBUG, FN_TYPE, NSCAT, JMUEO, MAXMSTRY
      PARAMETER(JMOUT = 6)
      COMMON / JMPARM /  PTJIM, YGAMMA, JMZMIN, JMRAD(264)
     &     ,PHAD, JMU2, JMV2, JMARRY( 6+MAXMS,0:NPSIMP )
     &     ,NLOST, TOTSCAT, ANOMOFF, JCMVAR, JMUEO
     &     ,JMPTYP(NPROC), JMBUG, FN_TYPE, MSFLAG, MAXMSTRY
      DOUBLE PRECISION JMPROC, JMVETO
      COMMON / JMEVNT/ JMPROC(NPROC)
     &,        JMVETO(2,13), NSCAT






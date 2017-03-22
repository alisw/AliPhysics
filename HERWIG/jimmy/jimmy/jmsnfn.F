      DOUBLE PRECISION FUNCTION JMSNFN(NDIM,V)
C ----------------------------------------------------------------------------
C  Purpose: Used when eikonalising the total hadronic xsec.
C ----------------------------------------------------------------------------

#include "herwig65.inc"
#include "jimmy.inc"


      INTEGER NDIM, N
      DOUBLE PRECISION XSECN, V(NDIM), B2MAX, B2, FAC1
      DOUBLE PRECISION CHI, JMAREA, JMFACT
      COMMON/BLOCK2/XSECN, N

      B2MAX = 5.D2
      FAC1  = B2MAX
      B2    = FAC1*V(1)

      CHI = JMAREA(B2)*XSECN*PHAD

      IF (JMBUG.GT.10) THEN
         WRITE(*,*) 'JMSNFN:XSECN,CHI',XSECN,CHI
      ENDIF

      IF (JMUEO.EQ.0) THEN
c     Standard (i.e. useless!) JIMMY. Generating QCD 2->2 events
c     totally democratically.
         JMSNFN = PIFAC/PHAD*FAC1*CHI**N/JMFACT(N)*DEXP(-CHI)
      ELSE
c     JIMMY being used to simulate the underlying event to a special
c     (low cross section) process.
         JMSNFN = N*PIFAC/PHAD*FAC1*CHI**N/JMFACT(N)*DEXP(-CHI)         
      ENDIF

      IF (JMBUG.GT.10) THEN
         WRITE(*,*) 'JMSNFN:N,JMFACT(N),JMSNFN',N,JMFACT(N),JMSNFN
      ENDIF

      RETURN
      END

#include "isajet/pilot.h"
C----------------------------------------------------------------------
      FUNCTION SUALFE(QS)
C----------------------------------------------------------------------
C
C     Returns the running EM coupling alpha_em(q**2)
C
C-----SEE BARGER/PHILLIPS, P. 202 ---------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sslun.inc"
      REAL SUALFE
      REAL PI,MB,SUM,A0,QD,QU,MS,ME,QS,MM,MD,MU,MTAU,MC
      DATA ME/.0005/,MM/.105/,MTAU/1.784/,MU/.01/,MD/.01/
      DATA MS/.5/,MC/1.6/,MB/5.0/,PI/3.1415926/
C
      SUM=0.
      QU=2./3.
      QD=-1./3.
      IF (QS.GT.4*ME**2) SUM=SUM+LOG(QS/4./ME**2)
      IF (QS.GT.4*MM**2) SUM=SUM+LOG(QS/4./MM**2)
      IF (QS.GT.4*MTAU**2) SUM=SUM+LOG(QS/4./MTAU**2)
      IF (QS.GT.4*MU**2) SUM=SUM+3*QU**2*LOG(QS/4./MU**2)
      IF (QS.GT.4*MD**2) SUM=SUM+3*QD**2*LOG(QS/4./MD**2)
      IF (QS.GT.4*MS**2) SUM=SUM+3*QD**2*LOG(QS/4./MS**2)
      IF (QS.GT.4*MC**2) SUM=SUM+3*QU**2*LOG(QS/4./MC**2)
      IF (QS.GT.4*MB**2) SUM=SUM+3*QD**2*LOG(QS/4./MB**2)
      A0=1./137.
      SUALFE=A0/(1.-A0/3./PI*SUM)
      RETURN
      END

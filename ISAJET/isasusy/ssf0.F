#include "isajet/pilot.h"
      COMPLEX*16 FUNCTION SSF0(XQSQ,XM1,XM2)
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      REAL XQSQ,XM1,XM2
      DOUBLE PRECISION QSQ,M1,M2,M1SQ,M2SQ,AQSQ,RE,XI,R
      DOUBLE PRECISION PI,T1,T2,BETA,XL,T
      DATA PI/3.14159265359D0/
      QSQ=XQSQ
      M1=XM1
      M2=XM2
      M1SQ=M1*M1
      M2SQ=M2*M2
      AQSQ=ABS(QSQ)
      IF(AQSQ.LT.1.D-6*(M1SQ+M2SQ)) THEN
        IF(ABS(M1-M2).LT.1.D-6*M1) THEN
          IF(M1SQ.LT.1.D-8) THEN
            RE=LOG(1.D-8)
          ELSE
            RE=LOG(M1SQ)+.5*(M2SQ/M1SQ-1.D0) -QSQ/(6.*M1SQ)
          ENDIF
        ELSE
          IF(M2SQ.LT.1.D-6*M1SQ) THEN
            RE=LOG(M1SQ)-1.D0
          ELSE IF(M1SQ.LT.1.D-6*M2SQ) THEN
            RE=LOG(M2SQ)-1.D0
          ELSEIF(M1SQ.GE.1.D-9*M2SQ) THEN
            RE=LOG(M1*M2)+(M1SQ+M2SQ)/(M1SQ-M2SQ)*LOG(M1/M2)-1.D0
     $      -QSQ*(.5*(M1SQ+M2SQ)-M1SQ*M2SQ*LOG(M2SQ/M1SQ)/(M2SQ-M1SQ))
     $      /(M1SQ-M2SQ)**2
          ENDIF
        ENDIF
        XI=0.D0
      ELSE
        IF(M1SQ.LT.1.D-6*AQSQ.OR.M2SQ.LT.1.D-6*AQSQ.OR.
     $  M1SQ+M2SQ.LT.1.D-5*AQSQ) THEN
          IF(M1SQ.LT.1.D-6*AQSQ) THEN
            R=M2SQ/QSQ
          ELSE
            R=M1SQ/QSQ
          ENDIF
          IF(ABS(R-1.D0).GT.1.D-6.AND.ABS(R).GT.1.D-6.AND.
     $    M1SQ+M2SQ.GT.1.D-5*AQSQ) THEN
            RE=LOG(R*QSQ)-2.D0+(1.D0-R)*LOG(ABS(1.D0-1.D0/R))
          ELSE
            RE=-2.D0+LOG(AQSQ)
          ENDIF
          IF(R.LT.1.D0.AND.R.GT.-1.D-10) THEN
            XI=-PI*(1.D0-R)
          ELSE
            XI=0.D0
          ENDIF
        ELSE
          T1=(M1-M2)*(M1-M2)
          T2=T1+4.0*M1*M2
          BETA=SQRT(ABS((1.D0-T1/QSQ)*(1.D0-T2/QSQ)))
          IF(QSQ.GT.T2.OR.QSQ.LT.T1) THEN
            XL=.5*BETA*LOG((QSQ*(1.D0+BETA)-M1SQ-M2SQ)/
     $      (QSQ*(1.D0-BETA)-M1SQ-M2SQ))
          ELSE
          T=M1SQ+M2SQ-QSQ
          IF(T.EQ.0.) T=1.D-10
            XL=BETA*ATAN(QSQ*BETA/T)
            IF(T.LT.0.D0) XL = XL + PI*BETA
          ENDIF
          RE=LOG(M1*M2)-(M1SQ-M2SQ)/QSQ*LOG(M2/M1)-2.D0+XL
          IF(QSQ.GT.T2) THEN
            XI=-PI*BETA
          ELSE
            XI=0.D0
          ENDIF
        ENDIF
      ENDIF
      SSF0=RE*(1.D0,0.D0)-XI*(0.D0,1.D0)
      RETURN
      END

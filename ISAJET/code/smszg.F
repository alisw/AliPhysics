#include "isajet/pilot.h"
      FUNCTION SMSZG(Q,QB,K,E,EB,AQ,BQ,AE,BE)
      IMPLICIT NONE
C
C     This does squared matrix element for q+qb -> Z+gamma
C          where Z-> e+eb
C     I have factored out 128*e^6*Q_q^2*|D_Z(z^2)| from
C       the expression. Also 1/12 from spin/color ave. is out.
C
      REAL Q(5),QB(5),K(5),E(5),EB(5),AQ,BQ,AE,BE,SMSZG
      REAL M1S,M2S,M12
      REAL EDQ,EBDK,EBDQ,EDK,QBDK,EDQB,EBDQB,QDK,QDQB
      EDQ=E(4)*Q(4)-E(1)*Q(1)-E(2)*Q(2)-E(3)*Q(3)
      EBDK=EB(4)*K(4)-EB(1)*K(1)-EB(2)*K(2)-EB(3)*K(3)
      EBDQ=EB(4)*Q(4)-EB(1)*Q(1)-EB(2)*Q(2)-EB(3)*Q(3)
      EDK=E(4)*K(4)-E(1)*K(1)-E(2)*K(2)-E(3)*K(3)
      QBDK=QB(4)*K(4)-QB(1)*K(1)-QB(2)*K(2)-QB(3)*K(3)
      EDQB=E(4)*QB(4)-E(1)*QB(1)-E(2)*QB(2)-E(3)*QB(3)
      EBDQB=EB(4)*QB(4)-EB(1)*QB(1)-EB(2)*QB(2)-EB(3)*QB(3)
      QDK=Q(4)*K(4)-Q(1)*K(1)-Q(2)*K(2)-Q(3)*K(3)
      QDQB=Q(4)*QB(4)-Q(1)*QB(1)-Q(2)*QB(2)-Q(3)*QB(3)
      M1S=(((AQ**2+BQ**2)*(AE**2+BE**2)-4*AQ*BQ*AE*BE)*EDQ*EBDK+
     $((AQ**2+BQ**2)*(AE**2+BE**2)+4*AQ*BQ*AE*BE)*EBDQ*EDK)/
     $4./QBDK
      M2S=(((AQ**2+BQ**2)*(AE**2+BE**2)+4*AQ*BQ*AE*BE)*EDQB*EBDK+
     $((AQ**2+BQ**2)*(AE**2+BE**2)-4*AQ*BQ*AE*BE)*EBDQB*EDK)/
     $4./QDK
      M12=(2*(AQ**2+BQ**2)*(AE**2+BE**2)*(EDQ*EBDQ*QBDK+EDQB*EBDQB*QDK)
     $+((AQ**2+BQ**2)*(AE**2+BE**2)-4*AQ*BQ*AE*BE)*(2*EDQ*EBDQB*QDQB+
     $EDQ*EBDK*QDQB-EDK*EBDQB*QDQB+EDQ*EBDQB*QDK-EDQ*EBDQB*QBDK)+
     $((AQ**2+BQ**2)*(AE**2+BE**2)+4*AQ*BQ*AE*BE)*(2*EDQB*EBDQ*QDQB-
     $EDQB*EBDK*QDQB+EDK*EBDQ*QDQB+EDQB*EBDQ*QDK-EDQB*EBDQ*QBDK))/
     $4./QBDK/QDK
      SMSZG=M1S+M2S+M12
      RETURN
      END

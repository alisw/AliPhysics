*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_UNIX)
      FUNCTION RNORTH (K)
C
C     Function used by NORRAN (CERN library entry V101) rewritten
C     for Alliant.                                     / Ch.Walck 871214
C
      DIMENSION C(45)
      EXTERNAL VNI, UNI
      DATA I1 /-71085056/, I2 /-25595858/
      DATA C / 0.9889430404, 0.9889430404, 0.9791515470, 0.9595685005
     +       , 0.9497770071, 0.9301939607, 0.9008194208, 0.8812363744
     +       , 0.8518618345, 0.8224872947, 0.7833212614, 0.7539466619
     +       , 0.7147806287, 0.6756145358, 0.6364485025, 0.5972824097
     +       , 0.5679078698, 0.5287418365, 0.4895757437, 0.4504097104
     +       , 0.4210351706, 0.3818690777, 0.3524945378, 0.3231199980
     +       , 0.2937454581, 0.2643709183, 0.2349963784, 0.2154133320
     +       , 0.1860387921, 0.1664557457, 0.1468726993, 0.1272896528
     +       , 0.1174981594, 0.0979151130, 0.0881236196, 0.0783321261
     +       , 0.0685405731, 0.0587490946, 0.0489575788, 0.0391660631
     +       , 0.0293745473, 0.0293745473, 0.0195830315, 0.0195830315
     +       , 0.0195830315 /
C
      IF ( K .GT. I1 ) GO TO 3
      S = UNI(0)
      T = UNI(0)
      B = AINT ( 7.0*(S+T) + 37.0*ABS(S-T) )
      X = UNI(0) - UNI(1)
      RNORTH = 0.0625 * ( X + SIGN (B,X) )
      RETURN
    3 IF ( K .GT. I2 ) GO TO 5
    4 RNORTH = 2.75 * VNI(0)
      J = 16.0 * ABS (RNORTH) + 1.0
      IF ( J - 14 ) 6, 6, 7
    6 P = (J+J-1) * .1497466E-2
      GO TO 8
    7 P = (89-J-J) * .698817E-3
    8 IF ( UNI(0) .GT. 79.78846 * (EXP (-0.5 * RNORTH * RNORTH)
     1       - C(J) - P * ( J - 16.0*ABS(RNORTH) ) ) ) GO TO 4
      RETURN
    5 V = VNI(0)
      IF ( V .EQ. 0.0 ) GO TO 5
      X = SQRT ( 7.5625 - 2.0*LOG(ABS(V)) )
      IF ( UNI(0)*X .GT. 2.75 ) GO TO 5
      RNORTH = SIGN ( X, V )
      RETURN
      END
#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:50  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE KMRSEB(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GL)
C***************************************************************C
C                                                               C
C     -----  VARIABLE QUARKS AND GLUONS AT SMALL X ----         C
C                                                               C
C     NEW VERSIONS !!!! JULY 1990                               C
C     "........................ " J. KWIECINSKI, A.D. MARTIN,   C
C     R.G. ROBERTS AND W.J. STIRLING PREPRINT DTP-90-46 )       C
C                                                               C
C  MODE 1 CORRESPONDS TO  HARRIMAN,                             C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  WITH SMALL X BEHAVIOUR DETERMINED FROM FIT  "HB FIT"         C
C                                                               C
C  MODE 2 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B0 FIT"         C
C                                                               C
C  MODE 3 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> X**-1/2 AS X--> 0 AT Q0**2    "B- FIT"         C
C                                                               C
C  MODE 4 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> X**-1/2 AS X--> 0 AT Q0**2    "B-(5) FIT"      C
C  I.E. WITH WEAK (R=5 GEV^-1) SHADOWING INCLUDED               C
C                                                               C
C  MODE 5 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 190 MEV  C
C  AND XG,XQ --> X**-1/2 AS X--> 0 AT Q0**2    "B-(2) FIT"      C
C  I.E. WITH STRONG  (R=2 GEV^-1) SHADOWING INCLUDED            C
C                                                               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C    THE FIRST NUMBER IN THE "HB" GRID IS  .03058               C
C    THE FIRST NUMBER IN THE "B0" GRID IS  .01727               C
C    THE FIRST NUMBER IN THE "B-" GRID IS  .01543               C
C    THE FIRST NUMBER IN THE "B-5"GRID IS  .0                   C
C    THE FIRST NUMBER IN THE "B-2"GRID IS  .0                   C
C                                                               C
C                                                               C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE LIGHT QUARK SEA I.E. UBAR(X)=DBAR(X)                C
C    = SEA/X FOR A PROTON.  IF IN DOUBT, CHECK THE              C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV)          C
C                                                               C
C                         -*-                                   C
C                                                               C
C     (THE RANGE OF APPLICABILITY IS CURRENTLY:                 C
C     10**-5 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO WJS@UK.AC.DUR.HEP             C
C                                                               C
C                                                               C
C***************************************************************C
#include "pdf/impdp.inc"
      IF(MODE.EQ.1) CALL STRH3B(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.2) CALL STRKB0O(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.3) CALL STRKBM(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.5) CALL SFKBMR2(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.5) CALL SFKBMR5(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      RETURN
      END

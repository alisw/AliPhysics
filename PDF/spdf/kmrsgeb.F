*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:51  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE KMRSGEB(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GL)
C***************************************************************C
C                                                               C
C     -----  VARIABLE LAMBDA AND GLUONS  ----                   C
C                                                               C
C     NEW VERSIONS !!!! OCTOBER 1990                            C
C     "........................ "  A.D. MARTIN,                 C
C     R.G. ROBERTS AND W.J. STIRLING PREPRINT DTP-90-?? )       C
C                                                               C
C  MODE 1 CORRESPONDS TO  KWIECINSKI,                           C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)                        C
C  WITH LAMBDA(4) = 190 MEV, ETAG = 5.10                        C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B0 FIT"         C
C                                                               C
C  MODE 2 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)                        C
C  WITH LAMBDA(4) = 135 MEV, ETAG = 4.65                        C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B135 FIT"       C
C                                                               C
C  MODE 3 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)                        C
C  WITH LAMBDA(4) = 160 MEV, ETAG = 4.25                        C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B160 FIT"       C
C                                                               C
C  MODE 4 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)                        C
C  WITH LAMBDA(4) = 200 MEV, ETAG = 5.65                        C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B200 FIT"       C
C                                                               C
C  MODE 5 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)                        C
C  WITH LAMBDA(4) = 235 MEV, ETAG = 5.20                        C
C  AND XG,XQ --> CONSTANT AS X--> 0 AT Q0**2   "B235 FIT"       C
C                                                               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C    THE FIRST NUMBER IN THE "B0"   GRID IS  .01727             C
C    THE FIRST NUMBER IN THE "B135" GRID IS  .01683             C
C    THE FIRST NUMBER IN THE "B160" GRID IS  .01663             C
C    THE FIRST NUMBER IN THE "B200" GRID IS  .01601             C
C    THE FIRST NUMBER IN THE "B235" GRID IS  .01571             C
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
      IF(MODE.EQ.1) CALL STRKB0O(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.2) CALL STRKB02(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.3) CALL STRKB03(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.4) CALL STRKB04(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.5) CALL STRKB05(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      RETURN
      END

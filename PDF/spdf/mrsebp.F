*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:54  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:25  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE MRSEBP(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GL)
C***************************************************************C
C                                                               C
C                                                               C
C         !!!!!!! NEW VERSIONS JANUARY 1989  !!!!!!!!           C
C                                                               C
C                                                               C
C  MODE 1 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (EMC FIT)    WITH LAMBDA= 100 MEV  C
C                                                               C
C  MODE 2  CORRESPONDS TO                                       C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 200 MEV  C
C                                                               C
C  (  SOFT GLUE :  X G(X,Q0) = A (1-X)**5 )                     C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE LIGHT QUARK SEA I.E. UBAR(X)=DBAR(X)=  ...          C
C    = SEA/X FOR A PROTON.  IF IN DOUBT, CHECK THE              C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV)          C
C                                                               C
C                         -*-                                   C
C                                                               C
C     (THE RANGE OF APPLICABILITY IS CURRENTLY:                 C
C     10**-4 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO SRG$T3@GEN                    C
C                                                               C
C                                                               C
C***************************************************************C
#include "pdf/impdp.inc"
      IF(MODE.EQ.1) CALL STRCEP(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.2) CALL STRCBP(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GL)
      RETURN
      END

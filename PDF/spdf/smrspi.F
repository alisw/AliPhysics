*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:30:53  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:48  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE SMRSPI(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C***************************************************************C
C                                                               C
C    PION DISTRIBUTION FUNCTIONS EXTRACTED FROM DRELL-YAN AND   C
C    PROMPT PHOTON PION BEAM DATA, USING (NUCLEAR TARGET        C
C    CORRECTED) HMRS(B) DISTRIBUTIONS FOR THE NUCLEON.          C
C    MAIN UNCERTAINTY IS IN SEA, HENCE FOLLOWING THREE SETS...  C
C    MODE=2 IS THE "BEST FIT" SET.                              C
C                                                               C
C    MODE 1 CORRESPONDS TO A 10% MOMENTUM SEA    (0.11213)      C
C    MODE 2 CORRESPONDS TO A 15% MOMENTUM SEA    (0.16119)      C
C    MODE 3 CORRESPONDS TO A 20% MOMENTUM SEA    (0.23785)      C
C                                                               C
C    (NUMBER IN BRACKETS IS LAST ENTRY IN FIRST ROW OF GRID)    C
C                                                               C
C                         -*-                                   C
C                                                               C
C    REFERENCE: A.D. MARTIN. R.G. ROBERTS. W.J. STIRLING        C
C    AND P.J. SUTTON, PREPRINT RAL-91-058 (1991).               C
C                         -*-                                   C
C                                                               C
C    FOR THESE SETS....                                         C
C    * LAMBDA(MSBAR,NF=4) = 190 MEV                             C
C    * UV=DV, STR=USEA=DSEA, BTM=0                              C
C                                                               C
C    NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION         C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE QUARK SEA / 6   I.E. UBAR(X)=DBAR(X)                C
C    = SEA/X FOR A PION. IF IN DOUBT, CHECK THE                 C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV           C
C                                                               C
C                         -*-                                   C
C                                                               C
C     THE RANGE OF APPLICABILITY IS CURRENTLY:                  C
C     10**-5 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO WJS@UK.AC.DUR.HEP             C
C                                                               C
C                                                               C
C***************************************************************C
#include "pdf/impdp.inc"
      IF(MODE.EQ.1) CALL SMRSP31(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.2) CALL SMRSP32(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.3) CALL SMRSP33(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      RETURN
      END

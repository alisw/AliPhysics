*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:52  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:24  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE MRSDSHD(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GL)
C***************************************************************C
C                                                               C
C                                                               C
C      MRS *** DIS SCHEME ** VERSIONS: S0',D0',D-' (FEB 1993)   C
C        NEW SET  MRS(H) ADDED    (DECEMBER 1993)               C
C                                                               C
C         MODE 1 IS THE  MRS(S0'-DIS) SET                       C
C         MODE 2 IS THE  MRS(D0'-DIS) SET                       C
C         MODE 3 IS THE  MRS(D-'-DIS) SET                       C
C         MODE 4 IS THE  MRS(H-DIS) SET                         C
C         ALL   SETS ARE FITTED TO THE RECENT NMC               C
C         AND CCFR  STRUCTURE FUNCTION DATA,                    C
C         AND HAVE (D0',D-',H)  UBAR < DBAR IN THE SEA.         C
C         BOTH THESE NEW SETS HAVE                              C
C         LAMBDA(MSbar,NF=4) = 230 MeV                          C
C                                                               C
C     THE REFERENCE IS: A.D. Martin, R.G. Roberts and           C
C     W.J. Stirling, RAL  preprint RAL-93-077 (1993)            C
C                                                               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C    THE FIRST NUMBER IN THE "1" GRID IS 0.01733                C
C    THE FIRST NUMBER IN THE "2" GRID IS 0.00930                C
C    THE FIRST NUMBER IN THE "3" GRID IS 0.00783                C
C    THE FIRST NUMBER IN THE "4" GRID IS 0.01508                C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC. IF IN DOUBT, CHECK THE  C
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
      IF(MODE.EQ.1) CALL SFDISS0(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     +GL)
      IF(MODE.EQ.2) CALL SFDISD0(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     +GL)
      IF(MODE.EQ.3) CALL SFDISDM(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     +GL)
      IF(MODE.EQ.4)CALL SFDISH(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GL)
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:53  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:25  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE MRSDSHP(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GL)
C***************************************************************C
C                                                               C
C                                                               C
C     NEW MRS VERSIONS: S0', D0', D-'  (NOVEMBER  1992)         C
C         MODE 1 IS THE 1990 KMRS(B0) SET;                      C
C         MODES 2-4 ARE NEW SETS FITTED TO THE RECENT NMC       C
C         AND CCFR ***FINAL***  STRUCTURE FUNCTION DATA.        C
C         DIFFERENCES BETWEEN THESE AND THE APRIL 1992          C
C         VERSIONS ARE SMALL! THE THREE NEW SETS HAVE           C
C         LAMBDA(MSbar,NF=4) = 230 MeV                          C
C                                                               C
C         MODE 5 is a new (October 1993) preliminary "fit"      C
C         to the recent H1 + ZEUS F2 data at small x,           C
C         otherwise similar to D0', D-'. The exponent of the    C
C         gluon and sea at small x for this new set is -0.3     C
C         (cf. -0.5 for D-')                                    C
C                                                               C
C     THE REFERENCE FOR SETS (2-4) IS:                          C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         Phys. Lett. B306 (1993) 145                           C
C                                                               C
C        MODE 1: KMRS(B0) (Lambda(4) = 190 MeV)                 C
C        MODE 2: MRS(S0') (updated B0, Lambda(4) = 230 MeV)     C
C        MODE 3: MRS(D0') (... but with ubar not= dbar)         C
C        MODE 4: MRS(D-') (updated B-, ubar not= dbar)          C
C        MODE 5: MRS(H) (a preliminary "HERA" fit, Oct. 1993)   C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C    THE FIRST NUMBER IN THE "1" GRID IS 0.01727                C
C    THE FIRST NUMBER IN THE "2" GRID IS 0.01617                C
C    THE FIRST NUMBER IN THE "3" GRID IS 0.00820                C
C    THE FIRST NUMBER IN THE "4" GRID IS 0.00678                C
C    THE FIRST NUMBER IN THE "5" GRID IS 0.01399                C
C                                                               C
C    NOTE THE EXTRA ARGUMENT IN THIS SUBROUTINE MRSEB,          C
C    TO ALLOW FOR THE POSSIBILITY OF A *** DIFFERENT ***        C
C    UBAR AND DBAR SEA!                                         C
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
      IF(MODE.EQ.1)CALL STRKB0(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GL)
      IF(MODE.EQ.2)CALL SFMSBS0(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     +GL)
      IF(MODE.EQ.3)CALL SFMSBD0(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     +GL)
      IF(MODE.EQ.4)CALL SFMSBDM(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     +GL)
      IF(MODE.EQ.5)CALL SFMSBH(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GL)
      RETURN
      END

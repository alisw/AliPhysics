*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:58  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      INTEGER FUNCTION VOMAS(IP)
      DIMENSION VMASS(21)
      EQUIVALENCE (IPC,APC)
      EQUIVALENCE(MPC,BPC)
C--                        ,  NOUGHT  ,   GAMMA  ,    NU    , ELECTRON ,
C--  ,    MU    ,   PI0    ,   PI+-   ,   K+-    ,    K0    ,  PROTON  ,
C--  , NEUTRON  ,  LAMBDA  ,  SIGMA+  ,  SIGMA0  ,  SIGMA-  ,   XI0    ,
C--  ,   XI-    , DEUTERON ,  TRITON  ,   HE3++  ,   HE4++  ,          =
      DATA  VMASS/          0.        ,.1E-10    ,.1E-9     ,.000511006,
     +  .105659 ,.13501    ,.13960    ,.4938     ,.4980     ,.938256   ,
     +  .939550 ,1.11540   ,1.18941   ,1.1923    ,1.19708   ,1.3143    ,
     +  1.3208  ,1.87549   ,2.80879   ,2.80824   ,3.72718/
      IPC=IP
      IF (IPC.LE.0) GO TO 15
      IF (IPC.GT.4096) GO TO 15
      IPC=MOD(IPC,64)
      IF (IPC.EQ.0) GO TO 15
      IF (IPC.GT.49) GO TO 15
      IF (IPC.LE.21) GO TO 20
      IF (IPC.GE.42) GO TO 20
C--                 RETURN UNCHANGED IF NOT AN INTEGER IN RANGE
   15 VOMAS=IP
      RETURN
   20 MPC=MOD(IPC,32)
      APC=32-IPC
      APC=SIGN(VMASS(MPC),APC)
      VOMAS=IPC
      RETURN
#if (defined(CERNLIB_CDC))&&(defined(CERNLIB_F4))
         ENTRY IOFMAS
#endif
#if !defined(CERNLIB_CDC)||!defined(CERNLIB_F4)
         ENTRY IOFMAS(IP)
#endif
C        ENTRY IOFMAS  (IP)
      IPC=IP
      MPC=IUCOMP(ABS(APC),VMASS,21)
      IF(APC.LT.0.)        MPC=MPC+32
      VOMAS=MPC
      RETURN
      END

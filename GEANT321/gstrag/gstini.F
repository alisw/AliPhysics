*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GSTINI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Initialization of tables for energy loss straggling in       *
C.    *   thin gas layers                                              *
C.    *                                                                *
C.    *   ==> Called by : GPHYSI                                       *
C.    *       Authors   : K. Lassila-Perini, I. Gavrilenko, P. Nevski  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcstra.inc"
C...
      DIMENSION GAMTT(21),FN(300),EN(300)
      DIMENSION ENOLD(300)
C...
      DATA GAMTT
     +  /     1.1 ,     1.2 ,     1.3 ,     1.5 ,     1.8 ,
     +        2.0 ,     2.5 ,     3.0 ,     4.0 ,     7.0 ,
     +       10.0 ,    20.0 ,    40.0 ,    70.0 ,   100.0 ,
     +      300.0 ,   600.0 ,  1000.0 ,  3000.0 , 10000.0 ,
     +    50000.0 /
C...
      DO 10 I=1,21
         GAMLOG(I) = LOG(GAMTT(I))
   10 CONTINUE
*
      JPHXS = LQ(JPHOT-1)
* NZ = number of elements in the cross-section bank
      NZ    = Q(JPHXS+1)
* NI = number of intervals
      NI    = Q(JPHXS+2+3*NZ)
      NFACT = 140/NI
      EMAX  = LOG(DCUTE*1E9)
      DO 30 NNP = 1,NI
         EM(NNP) = LOG(1E3*Q(JPHXS+2+3*NZ+1+(NNP-1)*5))
         IF (NNP.NE.NI) THEN
            EM2 = LOG(1E3*Q(JPHXS+2+3*NZ+1+NNP*5))
         ELSE
            EM2 = EMAX
         ENDIF
         EESTE = (EM2-EM(NNP))/NFACT
         DO 20 NPNP=1,NFACT
            EMM((NNP-1)*NFACT+NPNP) = EM(NNP) + (NPNP-1)*EESTE
            IF (NPNP.EQ.1.AND.NNP.NE.1)
     +        EMM((NNP-1)*NFACT+NPNP) = EM(NNP) + .3*EESTE
            IF (EMM((NNP-1)*NFACT+NPNP).GT.EMAX) THEN
               NP = (NNP-1)*NFACT+NPNP -1
               GOTO 40
            ENDIF
   20    CONTINUE
   30 CONTINUE
   40 EPS   = .001
*
      CALL GSTXIN
*
*
* Create and fill the banks. The number of intervals is decreased
* in GSTTAB. It is supposed that the EN vector is the same for
* every GAM, if not en error message is printed.
*
      DO 80 I = 1,21
         GAM = GAMTT(I)
         CALL GSTTAB(GAM,NT,EN,FN)
         IF (I.EQ.1) THEN
            CALL MZBOOK(IXCONS,JTSTCO,JTSTRA,-1,'STCO',0,0,NT*21,3,-1)
            CALL MZBOOK(IXCONS,JTSTEN,JTSTRA,-2,'STEN',0,0,NT,3,-1)
            DO 50 II=1,NT
               ENOLD(II) = EN(II)
   50       CONTINUE
         ELSE
            DO 60 III=1,NT
               IF (ENOLD(III).NE.EN(III)) PRINT *,'ERROR in GSTINI     '
     +         //'              (straggling): Energy vector EN is '
     +         //'changing in GSTTAB'
   60       CONTINUE
         ENDIF
         DO 70 J=1,NT
            Q(JTSTCO+(J-1)*21+I) = FN(NT+1-J)
   70    CONTINUE
   80 CONTINUE
*
* ***
      DO 90 J=1,NT
         Q(JTSTEN+J) = EN(NT+1-J)
   90 CONTINUE
*
      NNWE = 0
      DO 110 J=1,NT
         DO 100 I=1,21
            ANL1 = Q(JTSTCO+(J-1)*21+I)
            IF(ANL1.LE.0.) THEN
               GO TO 120
            ENDIF
            Q(JTSTCO+(J-1)*21+I) = LOG(ANL1)
  100    CONTINUE
         NNWE = NNWE+1
  110 CONTINUE
  120 Q(JTSTRA+1) = NNWE
*
      END

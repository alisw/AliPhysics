*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:34  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GSYNGE (GAMMA,AM,HKGAUS,STEP,NTR,ETR,STR)
C
C.    ******************************************************************
C.    *                                                                *
C.    *  Generation of synchrotron radiation photons:                  *
C.    *     o sample the number of generated photons from              *
C.    *       Poisson distribution                                     *
C.    *     o sample the energy of each photon from the                *
C.    *       tables below (these numbers can be computed by           *
C.    *       the routine GSYNIT)                                      *
C.    *       YT(I)  = photon energy/ critical photon energy (e/ec)    *
C.    *       FTI(I) = int from YT(I) to inf (20) of dN**2/dxde        *
C.    *                where dN**2/dxde = int from e/ec to inf of      *
C.    *                                   K_5/3(x) dx                  *
C.    *                                                                *
C.    *   GAMMA      - Lorentz factor,                                 *
C.    *   H          - transverse component of the magnetic field      *
C.    *                to particle direction,                          *
C.    *   STEP       - current step in magnetic field,                 *
C.    *   NTR        - number of the synchrotron photons,              *
C.    *   ETR(k)     - photon energy,                                  *
C.    *   STR(k)     - photon coordinate.                              *
C.    *                                                                *
C.    *   ==> Called by : GSYNC                                        *
C.    *       Authors   : Igor Gavrilenko                              *
C.    *                                                                *
C.    ******************************************************************
C
C
#include "geant321/gcunit.inc"
#if defined(CERNLIB_NEVER)
      COMMON /SYNTAB / YT(54),FT(53),FTI(53)
#endif
C...
      DIMENSION ETR(100),STR(100),RNDM(200)
#if !defined(CERNLIB_NEVER)
      DIMENSION YT(54),FTI(53)
#endif
C...
      DATA FTI/0.9827082157, 0.9756910801, 0.9708271027, 0.9669673443,
     +        0.9637124538, 0.9608694911, 0.9583284855, 0.9560201764,
     +        0.9538975954, 0.9519274235, 0.9368366003, 0.9263846874,
     +        0.9180962443, 0.9111109376, 0.9050130844, 0.8995658755,
     +        0.8946199417, 0.8900741935, 0.8858568072, 0.8536163568,
     +        0.8313701749, 0.8137849569, 0.7990072370, 0.7861415148,
     +        0.7746775150, 0.7642934322, 0.7547715902, 0.7459571362,
     +        0.6792119741, 0.6339961290, 0.5988188982, 0.5696868896,
     +        0.5446703434, 0.5226698518, 0.5029919147, 0.4851673245,
     +        0.4688622653, 0.3516088724, 0.2802021503, 0.2298262119,
     +        0.1918447465, 0.1620831937, 0.1381656080, 0.1185975820,
     +        0.1023733467,     0.8878208697E-01, 0.2343492396E-01,
     +        0.7047536317E-02, 0.2237016102E-02, 0.7283322047E-03,
     +        0.2375631739E-03, 0.7417966117E-04, 0.1891316060E-04/
 
      DATA YT/1E-5,2E-5,3E-5,4E-5,5E-5,6E-5,7E-5,8E-5,9E-5,
     +        1E-4,2E-4,3E-4,4E-4,5E-4,6E-4,7E-4,8E-4,9E-4,
     +        1E-3,2E-3,3E-3,4E-3,5E-3,6E-3,7E-3,8E-3,9E-3,
     +        1E-2,2E-2,3E-2,4E-2,5E-2,6E-2,7E-2,8E-2,9E-2,
     +        .1,.2,.3,.4,.5,.6,.7,.8,.9,1.,2.,3.,4.,5.,6.,7.,8.,9./
C...
      NTR         = 0
      H = 0.1*HKGAUS
      IF(H.LE.0.) GO TO 999
C...
      RAD         = 333.33*AM*GAMMA/H
C...
      EC          = 2.96E-8*GAMMA**3/RAD
C...
      AN          = .01053*GAMMA/RAD*STEP
      CALL GPOISS(AN,N,1)
      N=MIN(N,100)
      IF(N.LE.0) GO TO 999
C...
      CALL GRNDM(RNDM,N*2)
      DO 40 I = 1,N
         K1 = 1
         K2 = 53
         R = RNDM(2*I-1)*(FTI(1)-FTI(53))+FTI(53)
C...
   10    K = (K1+K2)/2
         IF(R.NE.FTI(K)) GO TO 20
         E = YT(K)*EC
         GO TO 30
C...
   20    IF(R.GT.FTI(K)) THEN
            K2 = K
         ELSEIF(R.LT.FTI(K)) THEN
            K1 = K
         ENDIF
         IF(K2-K1.GT.1) GO TO 10
C...
         IF(YT(K2).LE.1.) THEN
            Y2 = LOG(YT(K2))
            Y1 = LOG(YT(K1))
            E = Y2+(R-FTI(K2))*(Y1-Y2)/(FTI(K1)-FTI(K2))
            E = EXP(E)*EC
         ELSE
            F2 = LOG(FTI(K2))
            F1 = LOG(FTI(K1))
            DR = LOG(R)
            E = YT(K2)+(DR-F2)*(YT(K1)-YT(K2))/(F1-F2)
            E = E*EC
         END IF
C...
   30    CONTINUE
C
         IF(E/EC.GT.20.) THEN
            WRITE(CHMAIL,'('' *** GSYNGE: Photons have e>20*ec, e= '',
     +      G10.4,'' keV, ec= '',G10.4,'' keV'')') E,EC
            CALL GMAIL(1,0)
            GO TO 40
         ENDIF
C...
         NTR = NTR+1
         ETR(NTR) = E*1.E-6
         STR(NTR) = STEP*RNDM(2*I)
   40 CONTINUE
C...
 999  END

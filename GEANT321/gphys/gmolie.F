*
* $Id$
*
* $Log$
* Revision 1.2  1996/05/02 08:17:11  ravndal
* correct spheric coordinate sampling
*
* Revision 1.1.1.1  1995/10/24 10:21:27  cernlib
* Geant
*
*
*    EX GMOL
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMOLIE(OMEGA,BETA2,DIN)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Computes MOLIERE multiple scattering for a particle      *
C.    *       with parameters VECT in common /GCTRAK/                  *
C.    *                                                                *
C.    *       This subroutine must be called with the correct values   *
C.    *       of the constants OMC & CHC which depend of the medium    *
C.    *                                                                *
C.    *       OMC AND CHC are computed at initialisation time (GMOLI)  *
C.    *       No lateral displacement of the particle with respect     *
C.    *       the incident direction is included.                      *
C.    *       No path length correction is included                    *
C.    *                                                                *
C.    *  Generation of multiple scattering according to                *
C.    *  Moliere theory corrected for finite angle scattering          *
C.    *                                                                *
C.    *  evolved from Cern library program MLR                         *
C.    *                                                                *
C.    *  OMEGA & CHIC are number of scatterings and critical angle     *
C.    *  of the medium for a given incident particle                   *
C.    *                                                                *
C.    *  COSTH and SINTH are the cosine and sine of the generated      *
C.    *  scattering angle between 0 and 180 degrees                    *
C.    *                                                                *
C.    * THRED(NA)=reduced angles of Moliere theory                     *
C.    *                                                                *
C.    * F0I(NA),F1I(NA),F2I(NA)= integrale of Moliere functions        *
C.    *                                                                *
C.    *  4 point continued fraction interpolation is used to invert    *
C.    *  the total distribution function integral                      *
C.    *                                                                *
C.    * XINT= argument value                                           *
C.    * ARG= argument values of the table (arg,val)                    *
C.    * VAL= function values of the table (arg,val)                    *
C.    * THRI= the resulting interpolated function value                *
C.    *                                                                *
C.    *    ==>Called by : GMOLIE                                       *
C.    *         Author M.S. Dixit NRCC Ottawa    *********             *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctrak.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
      PARAMETER (ENEPER = 2.7182818)
      DIMENSION DIN(3),RNDM(3)
      DIMENSION TINT(40),ARG(4),VAL(4),THRED(40),F0I(40),F1I(40),F2I(40)
      DATA THRED/
     +   0.00, 0.10, 0.20, 0.30
     +,  0.40, 0.50, 0.60, 0.70
     +,  0.80, 0.90, 1.00, 1.10
     +,  1.20, 1.30, 1.40, 1.50
     +,  1.60, 1.70, 1.80, 1.90
     +,  2.00, 2.20, 2.40, 2.60
     +,  2.80, 3.00, 3.20, 3.40
     +,  3.60, 3.80, 4.00, 5.00
     +,  6.00, 7.00, 8.00, 9.00
     +, 10.00,11.00,12.00,13.00/
      DATA F0I/
     +  0.000000E+00 ,0.995016E-02 ,0.392106E-01 ,0.860688E-01
     + ,0.147856E+00 ,0.221199E+00 ,0.302324E+00 ,0.387374E+00
     + ,0.472708E+00 ,0.555142E+00 ,0.632121E+00 ,0.701803E+00
     + ,0.763072E+00 ,0.815480E+00 ,0.859142E+00 ,0.894601E+00
     + ,0.922695E+00 ,0.944424E+00 ,0.960836E+00 ,0.972948E+00
     + ,0.981684E+00 ,0.992093E+00 ,0.996849E+00 ,0.998841E+00
     + ,0.999606E+00 ,0.999877E+00 ,0.999964E+00 ,0.999990E+00
     + ,0.999998E+00 ,0.999999E+00 ,0.100000E+01 ,0.100000E+01
     + ,0.100000E+01 ,0.100000E+01 ,0.100000E+01 ,0.100000E+01
     + ,1.,1.,1.,1./
      DATA F1I/
     +  0.000000E+00,0.414985E-02,0.154894E-01,0.310312E-01
     + ,0.464438E-01,0.569008E-01,0.580763E-01,0.468264E-01
     + ,0.217924E-01,-0.163419E-01,-0.651205E-01,-0.120503E+00
     + ,-0.178272E+00,-0.233580E+00,-0.282442E+00,-0.321901E+00
     + ,-0.350115E+00,-0.366534E+00,-0.371831E+00,-0.367378E+00
     + ,-0.354994E+00,-0.314803E+00,-0.266539E+00,-0.220551E+00
     + ,-0.181546E+00,-0.150427E+00,-0.126404E+00,-0.107830E+00
     + ,-0.933106E-01,-0.817375E-01,-0.723389E-01,-0.436650E-01
     + ,-0.294700E-01,-0.212940E-01,-0.161406E-01,-0.126604E-01
     + ,-0.102042E-01,-0.840465E-02,-0.704261E-02,-0.598886E-02/
      DATA F2I/
     +  0.0,0.121500E-01,0.454999E-01,0.913000E-01
     + ,0.137300E+00,0.171400E+00,0.183900E+00,0.170300E+00
     + ,0.132200E+00,0.763000E-01,0.126500E-01,-0.473500E-01
     + ,-0.936000E-01,-0.119750E+00,-0.123450E+00,-0.106300E+00
     + ,-0.732800E-01,-0.312400E-01,0.128450E-01,0.528800E-01
     + ,0.844100E-01,0.114710E+00,0.106200E+00,0.765830E-01
     + ,0.435800E-01,0.173950E-01,0.695001E-03,-0.809500E-02
     + ,-0.117355E-01,-0.125449E-01,-0.120280E-01,-0.686530E-02
     + ,-0.385275E-02,-0.231115E-02,-0.147056E-02,-0.982480E-03
     + ,-0.682440E-03,-0.489715E-03,-0.361190E-03,-0.272582E-03/
*
*     ------------------------------------------------------------------
*
*
* *** Compute Theta angle from Moliere distribution
*
      CHIC  = CHCMOL*SQRT(STMULS)/(GETOT*BETA2)
      COSTH=1.
      SINTH=0.
      IF(OMEGA.LE.ENEPER)GO TO 90
      CNST=LOG(OMEGA)
      B=5.
      DO 10 L=1,10
         IF(ABS(B).LT.1.E-10)THEN
            B=1.E-10
         ENDIF
         DB=-(B-LOG(ABS(B))-CNST)/(1.-1./B)
         B=B+DB
         IF(ABS(DB).LE.0.0001)GO TO 20
   10 CONTINUE
      GO TO 90
   20 CONTINUE
      IF(B.LE.0.)GO TO 90
      BINV = 1./B
      TINT(1) = 0.
      DO 30 JA=2,4
         TINT(JA)=F0I(JA)+(F1I(JA)+F2I(JA)*BINV)*BINV
   30 CONTINUE
      NMAX = 4
   40 CONTINUE
      CALL GRNDM(RNDM,3)
      XINT=RNDM(2)
      DO 50 NA=3,40
         IF(NA.GT.NMAX) THEN
            TINT(NA)=F0I(NA)+(F1I(NA)+F2I(NA)*BINV)*BINV
            NMAX=NA
         ENDIF
         IF(XINT.LE.TINT(NA-1)) GO TO 60
   50 CONTINUE
      IF(XINT.LE.TINT(40)) THEN
         NA=40
         GOTO 60
      ELSE
         TMP=1.-(1.-B*(1.-XINT))**5
         IF(TMP.LE.0.)GO TO 40
         THRI=5./TMP
         GO TO 80
      ENDIF
   60 CONTINUE
      NA = MAX(NA-1,3)
      NA3 = NA-3
      DO 70 M=1,4
         NA3M=NA3+M
         ARG(M)=TINT(NA3M)
         VAL(M)=THRED(NA3M)**2
   70 CONTINUE
      F=THRED(NA)*.02
      CALL GMOL4(THRI,XINT,VAL,ARG,F,IER)
   80 CONTINUE
      TH=CHIC*SQRT(ABS(B*THRI))
      IF(TH.GT.PI)GO TO 40
      SINTH=SIN(TH)
      TEST=TH*(RNDM(3))**2
      IF(TEST.GT.SINTH)GO TO 40
      COSTH=COS(TH)
      GOTO 100
   90 CONTINUE
      CALL GRNDM(RNDM,1)
*
* *** Calculate SINE and COSINE of a random angle between 0 and 360 deg
*
  100 PHI = RNDM(1)*TWOPI
*
      DIN(1) = SINTH*COS(PHI)
      DIN(2) = SINTH*SIN(PHI)
      DIN(3) = COSTH
*
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:37  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_ASHO)
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GIASHO
#include "geant321/gcbank.inc"
#include "geant321/gcmate.inc"
      DIMENSION E0ELM(100),NSELM(100)
      DIMENSION ZSELM(7,100),ESELM(7,100)
      DATA E0ELM /
     + 0.0204,0.0385,3*0.,0.0738,0.0978,0.1157,0.1248,0.1338,
     + 3*0.,0.1745,0.1791,2*0.,0.1816,12*0.,0.3018,0.2806,
     + 0.2906,2*0.,0.3408,12*0.,0.4823,4*0.,0.5088,46*0. /
      DATA NSELM /
     + 1,1,3*0,2,2,2,2,2,3*0,3,3,2*0,3,12*0,4,4,4,2*0,4,12*0,
     + 5,4*0,5,46*0/
      DATA ZSELM /
     + 1.,0.,0.,0.,0.,0.,0.,2.,0.,0.,0.,0.,0.,0.,21*0.,
     + 4.,2.,0.,0.,0.,0.,0.,5.,2.,0.,0.,0.,0.,0.,
     + 6.,2.,0.,0.,0.,0.,0.,7.,2.,0.,0.,0.,0.,0.,
     + 8.,2.,0.,0.,0.,0.,0.,21*0.,
     + 4.,8.,2.,0.,0.,0.,0.,5.,8.,2.,0.,0.,0.,0.,14*0.,
     + 8.,8.,2.,0.,0.,0.,0.,84*0.,
     + 3.,18.,8.,2.,0.,0.,0.,4.,18.,8.,2.,0.,0.,0.,
     + 5.,18.,8.,2.,0.,0.,0.,14*0.,
     + 8.,18.,8.,2.,0.,0.,0.,84*0.,
     + 5.,18.,16.,8.,2.,0.,0.,28*0.,
     + 8.,18.,18.,8.,2.,0.,0.,322*0./
*23456789_123456789_123456789_123456789_123456789_123456789_123456789_12
      DATA ESELM / 0.01360,0.,0.,0.,0.,0.,0.,0.02459,0.,0.,0.,0.,0.,0.,
     +21*0., 0.01367,0.288,0.,0.,0.,0.,0.,0.01662,0.4030,0.,0.,0.,0.,
     +0., 0.01742,0.5380,0.,0.,0.,0.,0.,0.02174,0.6940,0.,0.,0.,0.,0.,
     +0.02643,0.8701,0.,0.,0.,0.,0.,21*0., .01047,.1147,1.844,0.,0.,0.,
     +0.,.01247,.1467,2.148,0.,0.,0.,0., 14*0.,0.01845,0.2666,3.206,0.,
     +0.,0.,0.,84*0., 0.00899,0.04480,1.169,10.37,0.,0.,0., 0.01063,
     +0.06190,1.274,11.11,0.,0.,0., 0.01291,0.07953,1.384,11.87,0.,0.,
     +0.,14*0., 0.01676,0.1412,1.750,14.33,0.,0.,0.,84*0., 0.00720,
     +0.04012,0.5682,3.908,27.93,0.,0.,28*0., 0.01466,0.1006,0.8097,
     +5.030,34.570,0.,0.,322*0./
*
*-----------------------------------------------------------------------
*
      DIMENSION ZSMED(50),ESMED(50)
C-----------------------------------------------------------------------
C
      JMA   = LQ(JMATE-NMAT)
      JMIXT = LQ(JMA-5)
      NCOMP = ABS(Q(JMA+11))
      AMED = 0.
      ZMED = 0.
      E0CAL = 0.   !Sum of Z(i)*log(I(i))
      E0MED = 0.
      NSMED = 0
 
      DO 20 I = 1,NCOMP
         IF(NCOMP.GT.1) THEN
            IZ = Q(JMIXT+NCOMP+I)+0.5
            WEIGHT = Q(JMIXT+2*NCOMP+I)/Q(JMIXT+I)
            E0MED = E0MED + Q(JMIXT+NCOMP+I)*WEIGHT*LOG(E0ELM(IZ))
            ZMED = ZMED + Q(JMIXT+NCOMP+I)*WEIGHT
            AMED = AMED + Q(JMIXT+I)*WEIGHT
         ELSE
            IZ = Z+0.5
            E0MED = E0ELM(IZ)
            ZMED = Z
            AMED = A
         ENDIF
         DO 10 J = 1,NSELM(IZ)
            NSMED = NSMED + 1
            IF(NCOMP.GT.1) THEN
               ZSMED(NSMED) = ZSELM(J,IZ)*WEIGHT
            ELSE
               ZSMED(NSMED) = ZSELM(J,IZ)
            ENDIF
            ESMED(NSMED) = ESELM(J,IZ)
            E0CAL = E0CAL + ZSMED(NSMED)*LOG(ESELM(J,IZ))
   10    CONTINUE
   20 CONTINUE
      IF (NCOMP.GT.1) E0MED = EXP(E0MED/ZMED)
      E0CAL = EXP(E0CAL/ZMED)
      ALFA = E0MED/E0CAL
C-----------------------------------------------------------------------
C      The following sets ZSMED and ESMED in the order of increase
C      of ESMED.
C-----------------------------------------------------------------------
      DO 40 I = 1,NSMED - 1
         IMIN = I
         EMIN = ESMED(I)
         DO 30 J = I + 1,NSMED
            IF (EMIN.LE.ESMED(J)) GOTO 30
            EMIN = ESMED(J)
            IMIN = J
   30    CONTINUE
         IF (I.EQ.IMIN) GOTO 40
         X = ESMED(I)
         Y = ZSMED(I)
         ESMED(I) = ESMED(IMIN)
         ZSMED(I) = ZSMED(IMIN)
         ESMED(IMIN) = X
         ZSMED(IMIN) = Y
   40 CONTINUE
C-----------------------------------------------------------------------
C      The following combines the first smallest oscillators whose
C      integer relative potentials are equal to 1.
C-----------------------------------------------------------------------
      ZMIN = 0.
      EMIN = 0.
      IMIN = 1
      DO 50 I = 1,NSMED
         J = ESMED(I)/ESMED(1) + 0.5
         IF (J.GT.1) GOTO 60
         ZMIN = ZMIN + ZSMED(I)
         EMIN = EMIN + ZSMED(I)*LOG(ESMED(I))
         IMIN = I
   50 CONTINUE
   60 ESMED(1) = EXP(EMIN/ZMIN)
      ZSMED(1)        = ZMIN
      DO 70 I = IMIN + 1,NSMED
         ZSMED(I - IMIN + 1) = ZSMED(I)
         ESMED(I - IMIN + 1) = ESMED(I)
   70 CONTINUE
      NSMED = NSMED - IMIN + 1
      JASHO = LQ(JMA-20)
*
* *** Store parameters of ASHO in material bank 20
      Q(JASHO+1) = NSMED
      Q(JASHO+2) = ZMED
      Q(JASHO+3) = AMED
      Q(JASHO+4) = ALFA
      Q(JASHO+5) = E0MED
      DO 80 KMED=1,NSMED
         Q(JASHO+5+KMED) = ZSMED(KMED)
         Q(JASHO+5+NSMED+KMED) = ESMED(KMED)
   80 CONTINUE
      NLEFT = 2*NSMED - 100
      CALL MZPUSH(IXCONS,JASHO,0,NLEFT,'I')
      END
#endif

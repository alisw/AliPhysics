*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 21/03/94  17.01.28  by  Rene Brun
*-- Author :
      SUBROUTINE GSCKOV(ITMED, NPCKOV, PPCKOV, ABSCO, EFFIC, RINDEX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Stores the tables for UV photon tracking in medium ITMED    *
C.    *    Please note that it is the user's responsability to         *
C.    *    provide all the coefficients:                               *
C.    *                                                                *
C.    *                                                                *
C.    *       ITMED       Tracking medium number                       *
C.    *       NPCKOV      Number of bins of each table                 *
C.    *       PPCKOV      Value of photon momentum (in GeV)            *
C.    *       ABSCO       Absorbtion coefficents                       *
C.    *                   dielectric: absorbtion length in cm          *
C.    *                   metals    : absorbtion fraction (0<=x<=1)    *
C.    *       EFFIC       Detection efficiency for UV photons          *
C.    *       RINDEX      Refraction index (if=0 metal)                *
C.    *                                                                *
C.    *       Called by : <USER>                                       *
C.    *                                                                *
C.    *       Authors: F.Carminati, R.Jones  ************              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcjloc.inc"
      DIMENSION PPCKOV(*), ABSCO(*), EFFIC(*), RINDEX(*)
C.
C.    ------------------------------------------------------------------
C.
      JTM = 0
      IF(ITMED.GT.0) THEN
         IF(JTMED.GT.0) THEN
            IF(IQ(JTMED-2).GE.ITMED) THEN
               JTM = LQ(JTMED-ITMED)
            ENDIF
         ENDIF
      ENDIF
      IF(JTM.LE.0) THEN
         WRITE(CHMAIL,10000)
10000   FORMAT(' **** GSCKOV Tracking medium ',I10,' not defined')
         CALL GMAIL(0,0)
         GO TO 999
      ENDIF
      IF(IQ(JTM-2).LT.3) THEN
         CALL MZPUSH(IXSTOR,JTM,3-IQ(JTM-2),0,' ')
      ENDIF
      CALL MZBOOK(IXCONS,JTCKOV,JTM,-3,'CKOV',4,4,NPCKOV+1,3,-1)
      CALL MZBOOK(IXCONS,JABSCO,JTCKOV,-1,'ABSC',0,0,NPCKOV,3,-1)
      CALL MZBOOK(IXCONS,JEFFIC,JTCKOV,-2,'EFFI',0,0,NPCKOV,3,-1)
      IF(RINDEX(1).GE.1.) THEN
         CALL MZBOOK(IXCONS,JINDEX,JTCKOV,-3,'RIND',0,0,NPCKOV,3,-1)
         CALL MZBOOK(IXCONS,JCURIN,JTCKOV,-4,'ABCU',0,0,NPCKOV,3,-1)
      ENDIF
      Q(JTCKOV+1)=NPCKOV
      DO 10 J=1, NPCKOV
         Q(JTCKOV+1+J) = PPCKOV(J)
         Q(JABSCO+  J) = MAX(1E-10,ABSCO(J))
         Q(JEFFIC+  J) = EFFIC(J)
         IF(RINDEX(1).GE.1.) THEN
            Q(JINDEX+ J) = RINDEX(J)
            IF(J.GT.1) THEN
               Q(JCURIN+J)=0.5*((1/RINDEX(J-1)**2)+(1./RINDEX(J)**2))
               Q(JCURIN+J)=Q(JCURIN+J-1)+(PPCKOV(J)-PPCKOV(J-1))*
     +                     Q(JCURIN+J)
            ELSE
               Q(JCURIN+J)=0.
            ENDIF
         ENDIF
   10 CONTINUE
*
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 06/07/94  17.30.05  by  S.Giani
*-- Author :
      SUBROUTINE GMXSEC(ILNMEC,SIGNN)
C
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcking.inc"
#include "geant321/gcnmec.inc"
C MICAP commons
#include "geant321/mmicap.inc"
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
      COMMON/MNUTRN/NAME,NAMEX,E,EOLD,NMED,MEDOLD,NREG,U,V,W,
     + UOLD,VOLD,WOLD,X,Y,ZZ,XOLD,YOLD,ZOLD,WATE,OLDWT,WTBC,
     + BLZNT,BLZON,AGE,OLDAGE,INEU,ENE(MAXNEU)
      INTEGER BLZNT
#include "geant321/mapoll.inc"
#include "geant321/mpoint.inc"
#include "geant321/mrecoi.inc"
#include "geant321/mmass.inc"
#include "geant321/mpstor.inc"
#include "geant321/cmagic.inc"
#include "geant321/mcreco.inc"
      DIMENSION IDCODE(NMECA-IBLOWN)
      DATA IDCODE /99,2,4,54,8,9,11,12,13,10,55,56,57,58,59,60,61,
     +            62,63,64,65,66/
C
C first check, if ZEBRA still in order
      IF(LD(LMAG1).NE.NMAGIC.OR.LD(LMAG2).NE.NMAGIC) THEN
         WRITE(6,*) ' CALOR: ZEBRA banks screwed up --> STOP'
         WRITE(IOUT,'('' MICAP: Magic number '',I12,'' not found: '',  '
     +   //'      2I12)') NMAGIC,LD(LMAG1),LD(LMAG2)
         STOP
      ENDIF
C Set the ID according to the process
      ID = IDCODE(ILNMEC-IBLOWN)
*
      E = GEKIN*1.E9
      EOLD = E
C Material number a la GEANT
*      NMED = NCEL
      NMED = NMAT
      NMEM=1
      CALL GTMED(NMED,IMED)
C get total cross-section
      CALL NSIGTA(E,NMED,TSIG,D,LD(LFP32),LD(LFP33))
C       THE PARAMETER (IIN) IS THE POINTER FOR ARRAYS DIMENSIONED BY
C       (NNUC) AND THE PARAMETER (IIM) IS THE POINTER FOR ARRAYS
C       DIMENSIONED BY (NMIX)
      SIGNN = 0.
      DO 10 IIM=1,NMIX
         IF(LD(LFP10+IIM-1).EQ.IMED) THEN
            IIN = LD(LFP16-1+IIM)
            CALL GMPLXS(D,LD, LD(LFP20),LD(LFP21),LD(LFP22),LD(LFP23),
     +      LD(LFP24), LD(LFP25),LD(LFP26),LD(LFP27),LD(LFP28),
     +      LD(LFP29), LD(LFP30), LD(LFP31),LD(LFP34),LD(LFP35),
     +      LD(LFP41),LD(LFP41+ NNUC), LD(LFP42),LD(LFP42+MEDIA),
     +      LD(LFP42+2*MEDIA),LD(LFP42+3* MEDIA), LD(LFP42+4*MEDIA),
     +      LD(LFP42+5*MEDIA),LD(LFP42+6*MEDIA), LD(LFP42+7*MEDIA),
     +      LD(LFP42+8*MEDIA),LD(LFP42+9*MEDIA), LD(LFP4 2+10*MEDIA),
     +      LD(LFP42+11*MEDIA),LD(LFP42+12*MEDIA), LD(LFP42+ 13*MEDIA),
     +      LD(LFP42+14*MEDIA),LD(LFP42+15*MEDIA), LD(LFP42+16* MEDIA),
     +      LD(LFP42+17*MEDIA),LD(LFP42+18*MEDIA), LD(LFP42+19* MEDIA),
     +      LD(LFP42+20*MEDIA),LD(LFP42+21*MEDIA), LD(LFP42+22* MEDIA),
     +      LD(LFP45),LD(LFP46),LD(LFP13), LD(LFP35+NQ*NNUC), D(LFP35+
     +      2*NQ*NNUC),IIN,IIM,ID,SIGNN1)
            SIGNN = SIGNN + SIGNN1*D(LFP12-1+IIM)
         ENDIF
   10 CONTINUE
C
      END

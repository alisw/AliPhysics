*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.47  by  S.Giani
*-- Author :
      SUBROUTINE GMICAP
C
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcking.inc"
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
C
C convert Z,A of recoil to CALOR particle code
C only p = 0, D = 7, T = 8, He3 = 9, alpha=10
* n = 13, p = 14, D = 45, T = 46, He3 = 49, alpha = 47
      DIMENSION NGPART(4,0:2)
      DATA ((NGPART(I,J),I=1,4),J=0,2)/13 ,-1 ,-1 ,
*                                      n
     +                               -1 , 14 , 45 ,
*                                         p    D
     +                               -1 , 46 , 49 ,
*                                         T    He3
     +                               -1 ,-1  , 47/
*                                              alpha
      SAVE
C first check, if ZEBRA still in order
      IF(LD(LMAG1).NE.NMAGIC.OR.LD(LMAG2).NE.NMAGIC) THEN
         WRITE(6,*) ' CALOR: ZEBRA banks screwed up --> STOP'
         WRITE(IOUT,'('' MICAP: Magic number '',I12,'' not found: '',  '
     +   //'      2I12)') NMAGIC,LD(LMAG1),LD(LMAG2)
         STOP
      ENDIF
C       THIS ROUTINE PERFORMS THE RANDOM WALK FOR ALL PARTICLES
   10 CONTINUE
C get material and particle information
*     U = UINC(1)
*     V = UINC(2)
*     W = UINC(3)
      U = VECT(4)
      V = VECT(5)
      W = VECT(6)
      X = 0.0
      Y = 0.0
      ZZ = 0.0
      BLZNT = 1
      WATE = 1.0
      AGE = 0.0
      NREG = 1
      WTBC = 1.0
C Energy MeV -> eV
*      E = EINC * 1.E6
      E = GEKIN*1.E9
C Material number a la GEANT
*      NMED = NCEL
      NMED = NMAT
      NMEM=1
C reset counter of heavy/charged and gamma bank
      NMEMR = 0
      NMEMG = 0
      INALB=0
      IET=0
      EOLD=E
      UOLD=U
      VOLD=V
      WOLD=W
      OLDWT=WATE
      XOLD=X
      YOLD=Y
      ZOLD=ZZ
      BLZON=BLZNT
      MEDOLD=NMED
      OLDAGE=AGE
      I=1
      CALL GTMED(NMED,IMED)
C get total cross-section
      CALL NSIGTA(E,NMED,TSIG,D,LD(LFP32),LD(LFP33))
C       DETERMINE WHICH ISOTOPE HAS BEEN HIT
      CALL ISOTPE(D,LD,LD(LFP10),D(LFP12),LD(LFP16),LD(LFP26),LD(LFP27),
     +            E,TSIG,IMED,IIN,IIM)
C       THE PARAMETER (IIN) IS THE POINTER FOR ARRAYS DIMENSIONED BY
C       (NNUC) AND THE PARAMETER (IIM) IS THE POINTER FOR ARRAYS
C       DIMENSIONED BY (NMIX)
      LD(LFP42+IMED-1)=LD(LFP42+IMED-1)+1
      INEU = 0
      NNEU = 0
      NHEVY = 0
      NGAMA = 0
      NPSTOR = 0
      CALL COLISN(D,LD,
     + LD(LFP20),LD(LFP21),LD(LFP22),LD(LFP23),LD(LFP24),
     + LD(LFP25),LD(LFP26),LD(LFP27),LD(LFP28),LD(LFP29),LD(LFP30),
     + LD(LFP31),D(LFP34),D(LFP35),LD(LFP41),LD(LFP41+NNUC),
     + LD(LFP42),LD(LFP42+MEDIA),LD(LFP42+2*MEDIA),LD(LFP42+3*MEDIA),
     + LD(LFP42+4*MEDIA),LD(LFP42+5*MEDIA),LD(LFP42+6*MEDIA),
     + LD(LFP42+7*MEDIA),LD(LFP42+8*MEDIA),LD(LFP42+9*MEDIA),
     + LD(LFP42+10*MEDIA),LD(LFP42+11*MEDIA),LD(LFP42+12*MEDIA),
     + LD(LFP42+13*MEDIA),LD(LFP42+14*MEDIA),LD(LFP42+15*MEDIA),
     + LD(LFP42+16*MEDIA),LD(LFP42+17*MEDIA),LD(LFP42+18*MEDIA),
     + LD(LFP42+19*MEDIA),LD(LFP42+20*MEDIA),LD(LFP42+21*MEDIA),
     + LD(LFP42+22*MEDIA),LD(LFP45),LD(LFP46),LD(LFP13),
     + LD(LFP35+NQ*NNUC),D(LFP35+2*NQ*NNUC),IIN,IIM)
      CALL BANKR(D,LD,5)
C -------- fill return arrays with generated particles ---------------
C first heavy/charged particles
   20 NPHETC = 0
      NRECOL = 0
      ERMED  = 0.0
      EETOT = 0.0
C -------- store  neutrons -------------------------------------
      INTCAL = 0
C
      ISTOP = 1
      JPA = LQ(JPART-13)
      AGEMNE = Q(JPA+7)
      NGKINE = 0
      DO 30  N=1,NNEU
         CALL GETPAR(IDNEU,N,IERR)
         IF(IERR.EQ.0) THEN
            NGKINE = NGKINE + 1
            TTKIN  = EP * 1.E-9
            PGEANT = SQRT(TTKIN*(TTKIN+2*AGEMNE))
            GKIN(1,NGKINE) = UP*PGEANT
            GKIN(2,NGKINE) = VP*PGEANT
            GKIN(3,NGKINE) = WP*PGEANT
            GKIN(4,NGKINE) = TTKIN + AGEMNE
            GKIN(5,NGKINE) = 13
            TOFD(NGKINE)   = AGEP * 1.E-9
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
*           NPHETC = NPHETC + 1
*           IF(NPHETC.GT.MXCP) NPHETC=MXCP
*           IPCAL(NPHETC) = 1
C kinetic energy in MeV
*           EKINET(NPHETC) = EP * 1.E-6
*           UCAL(NPHETC,1) = UP
*           UCAL(NPHETC,2) = VP
*           UCAL(NPHETC,3) = WP
*           CALTIM(NPHETC) = AGEP
         ENDIF
   30 CONTINUE
C -------- store heavy recoil products ------------------------
      DO 40  N=1,NHEVY
         CALL GETPAR(IDHEVY,N,IERR)
         IF(IERR.EQ.0) THEN
C check particle type
            MA = NINT(AMP)
            MZ = NINT(ZMP)
            IF(MA.LE.4.AND.MZ.LE.2) THEN
               IF(NGPART(MA,MZ).EQ.-1) GOTO 40
            ELSE
C get heavy recoil nucleus
               NRECOL = NRECOL + 1
               AMED(NRECOL) = AMP
               ZMED(NRECOL) = ZMP
               ERMED = ERMED + EP * 1.E-9
               GOTO 40
            ENDIF
C store particle type
            NGKINE = NGKINE + 1
            JPA = LQ(JPART-NGPART(MA,MZ))
            AGEMAS = Q(JPA+7)
            TTKIN  = EP * 1.E-9
            PGEANT = SQRT(TTKIN*(TTKIN+2*AGEMAS))
            GKIN(1,NGKINE) = UP*PGEANT
            GKIN(2,NGKINE) = VP*PGEANT
            GKIN(3,NGKINE) = WP*PGEANT
            GKIN(4,NGKINE) = TTKIN + AGEMAS
            GKIN(5,NGKINE) = NGPART(MA,MZ)
            TOFD(NGKINE) = AGEP * 1.E-9
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
*           NPHETC = NPHETC + 1
*           IF(NPHETC.GT.MXCP) NPHETC=MXCP
*           IPCAL(NPHETC) = NPART(MA,MZ)
C kinetic energy in MeV
*           EKINET(NPHETC) = EP * 1.E-6
*           UCAL(NPHETC,1) = UP
*           UCAL(NPHETC,2) = VP
*           UCAL(NPHETC,3) = WP
*           CALTIM(NPHETC) = AGEP
         ENDIF
   40 CONTINUE
* Number of produced particles (may be > MXGKIN)
      NNEHEG = NGKINE + NGAMA
C
C----------- get generated gammas --------------------
      NS   = 0
      NREM = 0
      DO 50  N=1,NGAMA
         IF (NS.GE.NGAMA) GO TO 60
         NS = NS + 1
         CALL GETPAR(IDGAMA,NS,IERR)
         IF(IERR.EQ.0) THEN
            IF (NNEHEG-NREM.GT.MXGKIN) THEN
               NREM  = NREM + 1
               UP1   = UP
               VP1   = VP
               WP1   = WP
               EP1   = EP
               AGEP1 = AGEP
               NS = NS + 1
*    Get the other gamma to be summed with the previous one
               CALL GETPAR(IDGAMA,NS,IERR)
               IF(IERR.EQ.0) THEN
                 UP  = (UP1*EP1+UP*EP)
                 VP  = (VP1*EP1+VP*EP)
                 WP  = (WP1*EP1+WP*EP)
*    Normalize the new direction cosines
                 WUP = SQRT(UP**2+VP**2+WP**2)
                 UP  = UP/WUP
                 VP  = VP/WUP
                 WP  = WP/WUP
                 EP  = EP1 + EP
                 AGEP = AGEP1 + AGEP
               ENDIF
            ENDIF
            NGKINE = NGKINE + 1
            PGEANT = EP * 1.E-9
            GKIN(1,NGKINE) = UP*PGEANT
            GKIN(2,NGKINE) = VP*PGEANT
            GKIN(3,NGKINE) = WP*PGEANT
            GKIN(4,NGKINE) = PGEANT
            GKIN(5,NGKINE) = 1
            TOFD(NGKINE) = AGEP * 1.E-9
            GPOS(1,NGKINE) = VECT(1)
            GPOS(2,NGKINE) = VECT(2)
            GPOS(3,NGKINE) = VECT(3)
*           NG = NG + 1
*           NPHETC = NPHETC + 1
*           IF(NPHETC.GT.MXCP) NPHETC=MXCP
*           IPCAL(NPHETC) = 11
*           EKINET(NPHETC) = EP*1.E-6
*           UCAL(NPHETC,1) = UP
*           UCAL(NPHETC,2) = VP
*           UCAL(NPHETC,3) = WP
*           CALTIM(NPHETC) = AGEP
C nucleus is in ground state !
            EXMED = 0.0
         ENDIF
   50 CONTINUE
* only one neutron generated -> the particle is the same
   60 IF (NGKINE.EQ.1.AND.GKIN(5,1).EQ.13) THEN
         NGKINE  = 0
         CALL GETPAR(IDNEU,1,IERR)
         VECT(4) = UP
         VECT(5) = VP
         VECT(6) = WP
         GEKIN   = EP * 1.E-9
         GETOT   = GEKIN + AGEMNE
         VECT(7) = SQRT(GEKIN*(GEKIN+2.*AGEMNE))
         TOFG    = TOFG + AGEP * 1.E-9
         ISTOP   = 0
      ENDIF
*
      IF (MTP        .EQ.         2) THEN
         INTCAL = 13
      ELSEIF (MTP    .EQ.        18) THEN
         IF (NHEVY.GT.0) INTCAL = 15
      ELSEIF (MTP    .LT.       100) THEN
         IF (NNEU .GT.0) INTCAL = 20
      ELSEIF (MTP    .EQ.       102) THEN
         IF (NGAMA.GT.0) INTCAL = 18
      ELSEIF (MTP    .GE.       100) THEN
         IF (NHEVY+NGAMA.GT.0) INTCAL = 16
      ENDIF
      IF(NNEU+NHEVY+NGAMA.GT.0.AND.INTCAL.EQ.0) INTCAL=12
      KCASE = NAMEC(INTCAL)
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:35  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE BLKDT4.FOR
*COPY BLKDT4
*
*=== blkdt4 ===========================================================*
*==                                                                    *
      SUBROUTINE FLKDT4
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*-----  Block Data without names
#include "geant321/hadflg.inc"
#include "geant321/redver.inc"
      COMMON / FKSLOP  / AMBMBB(75)
      DIMENSION IKII1(17), IEII1(17), IRII1(17)
      DIMENSION AMBMB1(75)
      DIMENSION IKCHX1(16), INNUR1(2,2,16)
*
*     datas     datas    datas      datas     datas
*------          *********
      DATA IKII1/ 0, 15, 41, 67, 82, 93, 110, 133, 148, 159, 172, 183,
     &           207, 224, 241, 252, 268 /
      DATA IEII1/ 0, 21, 46, 71, 92, 109, 126, 143, 160, 173, 186, 199,
     &           220, 241, 262, 279, 296 /
      DATA IRII1/ 0,315, 965, 1615, 1930, 2117, 2406, 2797, 3052, 3195,
     &           3364, 3507, 4011, 4368, 4725, 4912, 5184/
 
*
*     Masses for the slope B(M) in GeV
*     Slope B(M) for an mesonic system
*     Slope B(M) for a baryonic system
*
      DATA AMBMB1/  0.8D0, 0.85D0,  0.9D0, 0.95D0, 1.D0,
     &     1.05D0,  1.1D0, 1.15D0,  1.2D0, 1.25D0,
     &      1.3D0,  1.35D0, 1.4D0,  1.45D0,  1.5D0,
     &     1.55D0,  1.6D0,  1.65D0, 1.7D0,   1.75D0,
     &      1.8D0,  1.85D0, 1.9D0,  1.95D0,  2.D0,
     &     15.6D0, 14.95D0, 14.3D0, 13.65D0, 13.D0,
     &    12.35D0, 11.7D0, 10.85D0, 10.D0,  9.15D0,
     &      8.3D0,  7.8D0,  7.3D0,  7.25D0,  7.2D0,
     &     6.95D0,  6.7D0,  6.6D0,  6.5D0,   6.3D0,
     &      6.1D0,  5.85D0, 5.6D0,  5.35D0,  5.1D0,
     &      15.D0,   15.D0, 15.D0,  15.D0,   15.D0, 15.D0, 15.D0,
     &     14.2D0,  13.4D0, 12.6D0,
     &     11.8D0, 11.2D0, 10.6D0,  9.8D0,    9.D0,
     &     8.25D0,  7.5D0, 6.25D0,  5.D0,    4.5D0, 5*4.D0 /
*
*  /Hadflg/
      DATA IELFL1, ICXFL1, IOLDH1 / 3*0 /
      DATA IKCHX1 / 0, 2, 2, 0, 0, 2, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0 /
      DATA INNUR1 / 13, 1, 2*0, 14, 8, 23, 1, 14, 1, 23, 8, 14, 8, 2*0,
     &              15, 1, 2*0, 15, 8, 24, 1, 16, 1, 25, 8, 16, 8, 2*0,
     &               1, 1, 2*0,  1, 8,  8, 1,  8, 8, 2*0,  2, 1,  9, 8,
     &               2, 8, 2*0,  9, 1, 2*0, 24, 8, 2*0, 25, 1, 2*0 /
*
      DO 10 J=1,17
        IKII(J) = IKII1(J)
        IEII(J) = IEII1(J)
        IRII(J) = IRII1(J)
  10  CONTINUE
*
      DO 20 J=1,75
         AMBMBB(J) = AMBMB1(J)
   20 CONTINUE
*
      DO 30 J=1,16
         IKCHXG(J) = IKCHX1(J)
   30 CONTINUE
      DO 60 J1=1,2
        DO 50 J2=1,2
          DO 40 J3=1,16
            INNURE(J1,J2,J3) = INNUR1(J1,J2,J3)
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
      IELFLG = IELFL1
      ICXFLG = ICXFL1
      IOLDHD = IOLDH1
*=                                               end*block.blkdt4      *
      END

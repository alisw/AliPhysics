*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*FCA :          17/02/99  12:00:00  by  Federico Carminati
*               Calculation of radiation length and interaction length
*               ignoring user input parameters
*-- Author :
      SUBROUTINE GSMATE(IMAT,NAMATE,A,Z,DENS,RADL,ABSL,UBUF,NWBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       Store material parameters                                *
C.    *                                                                *
C.    *                                                                *
C.    *             The Material data structure JMATE                  *
C.    *             ---------------------------------                  *
C.    *                                                                *
C.    *                                         | JMATE                *
C.    *    NMATE           IMATE                v                      *
C.    *     ......................................                     *
C.    *     |              | |                 | |                     *
C.    *     ......................................                     *
C.    *                     |                                          *
C.    *                     | JMA                                      *
C.    *                     v                                          *
C.    *                    .....................                       *
C.    *                    | 1 |               |                       *
C.    *                    .....               |                       *
C.    *                    | 2 |  Material     |                       *
C.    *                    |...|               |                       *
C.    *                    | 3 |   Name        |                       *
C.    *                    |...|               |                       *
C.    *                    | 4 |               |                       *
C.    *                    |...|               |                       *
C.    *                    | 5 |               |                       *
C.    *                    .....................                       *
C.    *                    | 6 |   A           |                       *
C.    *                    |...|...............|                       *
C.    *                    | 7 |   Z           |                       *
C.    *                    |...|...............|                       *
C.    *                    | 8 |   Density     |                       *
C.    *                    |...|...............|                       *
C.    *                    | 9 |   RADL        |                       *
C.    *                    |...|...............|                       *
C.    *                    | 10|   ABSL        |                       *
C.    *                    |...|...............|                       *
C.    *                    | 11|   NMIXT       |                       *
C.    *                    |...|...............|                       *
C.    *                    |   |               |                       *
C.    *                    .....................                       *
C.    *                                                                *
C.    * JMA = LQ(JMATE-IMATE) pointer to material IMATE                *
C.    *                                                                *
C.    *   When  the subroutine  GPHYSI is  called at  initialisation   *
C.    * time  the following  banks  are  created for  each  material   *
C.    * (tabulation of energy loss and cross-section).                 *
C.    *                                               | JMATE          *
C.    * NMATE                     IMATE               v                *
C.    * ................................................               *
C.    * |                        | |                 | |               *
C.    * ................................................               *
C.    *                           | JMA = LQ(JMATE-IMATE)              *
C.    *                           v                              11    *
C.    * ............................................................   *
C.    * |  13 12 11 10 9 8 7 6 5 4 3 2 1 | |  Material parameters  |   *
C.    * ............................................................   *
C.    *    |  |  |  |  | | | | | | | | |                               *
C.    *    |  |  |  |  | | | | | | | | v  JMAEL = LQ(JMA-1)            *
C.    *    |  |  |  |  | | | | | | | |                         270     *
C.    *    |  |  |  |  | | | | | | | |................................ *
C.    *    |  |  |  |  | | | | | | | ||Energy loss for electron/positro*
C.    *    |  |  |  |  | | | | | | | |............................     *
C.    *    |  |  |  |  | | | | | | | v  JMAMU = LQ(JMA-2)       90     *
C.    *    |  |  |  |  | | | | | | |..............................     *
C.    *    |  |  |  |  | | | | | | ||Energy loss for muons       |     *
C.    *    |  |  |  |  | | | | | | |..............................     *
C.    *    |  |  |  |  | | | | | | v  JMAAL = LQ(JMA-3)         90     *
C.    *    |  |  |  |  | | | | | |................................     *
C.    *    |  |  |  |  | | | | | ||Energy loss for other particles|    *
C.    *    |  |  |  |  | | | | | |................................     *
C.    *    |  |  |  |  | | | | | v  JPROB = LQ(JMA-4)           30     *
C.    *    |  |  |  |  | | | | |..................................     *
C.    *    |  |  |  |  | | | | ||Some material constants         |     *
C.    *    |  |  |  |  | | | | |..................................     *
C.    *    |  |  |  |  | | | | v  JMIXT = LQ(JMA-5)             11     *
C.    *    |  |  |  |  | | | |....................................     *
C.    *    |  |  |  |  | | | ||Mixture or compound parameters    |     *
C.    *    |  |  |  |  | | | |....................................     *
C.    *    |  |  |  |  | | | v  JPHOT = LQ(JMA-6) and JMUNU     90     *
C.    *    |  |  |  |  | | |......................................     *
C.    *    |  |  |  |  | | ||Photo-effect cross-section          |     *
C.    *    |  |  |  |  | | |......................................     *
C.    *    |  |  |  |  | | v  JANNI = LQ(JMA-7)                 90     *
C.    *    |  |  |  |  | |........................................     *
C.    *    |  |  |  |  | ||Positron annihilation cross-section   |     *
C.    *    |  |  |  |  | |........................................     *
C.    *    |  |  |  |  | V  JCOMP = LQ(JMA-8)                   90     *
C.    *    |  |  |  |  |..........................................     *
C.    *    |  |  |  |  ||Compton scattering cross-section        |     *
C.    *    |  |  |  |  |..........................................     *
C.    *    |  |  |  |  V  JBREM = LQ(JMA-9)                     90     *
C.    *    |  |  |  | ............................................     *
C.    *    |  |  |  | |Bremsstrahlung cross-section              |     *
C.    *    |  |  |  | ............................................     *
C.    *    |  |  |  V  JPAIR = LQ(JMA-10)                       90     *
C.    *    |  |  | ...............................................     *
C.    *    |  |  | |Pair production cross-section                |     *
C.    *    |  |  | ...............................................     *
C.    *    |  |  V  JDRAY = LQ(JMA-11)                         210     *
C.    *    |  | ..................................................     *
C.    *    |  | |Moller and Bhabha cross-sections                |     *
C.    *    |  | ..................................................     *
C.    *    |  V  JPFIS = LQ(JMA-12)                             90     *
C.    *    | .....................................................     *
C.    *    | |Photo fission cross section                        |     *
C.    *    | .....................................................     *
C.    *    V  JRAYL = LQ(JMA-13)                                62     *
C.    *   ........................................................     *
C.    *   |Rayleigh scattering cross section and atomic form fact|     *
C.    *   ........................................................     *
C.    * V  JMUNU = LQ(JMA-14)                                  90      *
C.    *   ........................................................     *
C.    * V  JRANG = LQ(JMA-15)                                 180      *
C.    * V........................................................      *
C.    *  |Stopping range for electrons/positrons                |      *
C.    *  ........................................................      *
C.    * V  JRANG = LQ(JMA-16)                                 180      *
C.    * V........................................................      *
C.    *  |Stopping range for muons / other particles            |      *
C.    *  ........................................................      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGEOM    ,<GXINT> GINC3              *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gcunit.inc"
#include "geant321/gconsp.inc"
      DIMENSION UBUF(1)
      CHARACTER*(*) NAMATE
      CHARACTER*20 NAME
      PARAMETER(ALR2AV=1.39621E-03, AL183=5.20948)
C.
C.    ------------------------------------------------------------------
C.
      IF(IMAT.LE.0)GO TO 99
      IF(JMATE.LE.0)THEN
         CALL MZBOOK(IXCONS,JMATE,JMATE,1,'MATE',NMATE,NMATE,0,3,0)
         IQ(JMATE-5)=0
      ENDIF
      IF(IMAT.GT.NMATE)THEN
         CALL MZPUSH(IXCONS,JMATE,IMAT-NMATE,0,'I')
         NMATE=IMAT
         JMA1=0
      ELSE
         JMA1=LQ(JMATE-IMAT)
         IF(JMA1.GT.0) THEN
            WRITE(CHMAIL,10000)
            CALL GMAIL(1,0)
            CALL GPMATE(IMAT)
            CALL MZDROP(IXCONS,LQ(JMATE-IMAT),' ')
         ENDIF
      ENDIF
      CALL MZBOOK(IXCONS,JMA,JMATE,-IMAT,'MATE',20,20,NWBUF+11,IOMATE,0)
C
      NAME=NAMATE
      NCH=LNBLNK(NAME)
      IF(NCH.GT.0)THEN
         IF(NAME(NCH:NCH).EQ.'$')NAME(NCH:NCH)=' '
      ENDIF
      CALL UCTOH(NAME,IQ(JMA+1),4,20)
C
      Q(JMA + 6) = A
      Q(JMA + 7) = Z
      Q(JMA + 8) = DENS
* Calculation with standard formulas
*     Q(JMA + 9) = RADL
*     Q(JMA + 10) = ABSL
      IF(A.GT.0.AND.Z.GT.0) THEN
         IF(RADL.LT.0) THEN
            Q(JMA+9)=-RADL
         ELSE
            ALZ=LOG(Z)/3
            Q(JMA + 9) = A/(ALR2AV*DENS*Z*(Z+GXSI(Z))*
     +        (AL183-LOG(Z)/3-GFCOUL(Z)))
         ENDIF
         IF(ABSL.LT.0) THEN
            Q(JMA+10)=-ABSL
         ELSE
            Q(JMA + 10) = A/(AVO*DENS*1E-3*GHSIGM(5.,8,A))
         ENDIF
      ELSE
         Q(JMA + 9) = BIG
         Q(JMA + 10) = BIG
      ENDIF
      Q(JMA + 11) = 1.
      IF(NWBUF.GT.0)CALL UCOPY(UBUF,Q(JMA+12),NWBUF)
C
      IF(JMA1.GT.0) THEN
         CALL GPMATE(-IMAT)
      ENDIF
C
10000 FORMAT(' *** GSMATE ***: Warning, material redefinition:')
  99  END
 

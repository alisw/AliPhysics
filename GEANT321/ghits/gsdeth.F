*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:12  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GSDETH(IUSET,IUDET,NH,NAMESH,NBITSH,ORIG,FACT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *          Handling sensitive DETector Hit parameters            *
C.    *          ------------------------------------------            *
C.    *                                                                *
C.    *   Defines hit parameters for detector IUDET of set IUSET.      *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user detector identifier                             *
C.    * NH        number of elements per hit                           *
C.    * NAMESH    the NH variable names for the hit elements           *
C.    * NBITSH    the NH bit numbers for packing the variable values   *
C.    * ORIG      The quantity packed in the structure JHITS for the   *
C.    *           Ith variable is a  positive integer with NBITSH(I)   *
C.    *           bits and such that                                   *
C.    * FACT      IVAR(I) = (VAR(I)+ORIG(I))*FACT(I)                   *
C.    *   The  routine  is  used at  initialisation  time  once  the   *
C.    * geometrical volumes  have been defined  to describe  the hit   *
C.    * elements and the way to do packing in memory and on tape.      *
C.    *                           EXAMPLE                              *
C.    *   Assume an electromagnetic calorimeter ECAL divided into 40   *
C.    * PHI  sections called  EPHI.   Each  EPHI division  is  again   *
C.    * divided along  the Z axis in  60 objects called  EZRI.  Each   *
C.    * EZRI  is finally  divided into  4 lead  glass blocks  called   *
C.    * BLOC.                                                          *
C.    * The geometrical  information to describe  one hit  will then   *
C.    * be:                                                            *
C.    *  - The EPHI section number (between 1 and 40)                  *
C.    *  - The EZRI division number (between 1 and 60)                 *
C.    *  - The BLOC number (1 to 4)                                    *
C.    * The variables we want to store for each hit are for example:   *
C.    *  - X     x position of the hit in the lead glass block         *
C.    *  - Y     y                                                     *
C.    *  - Z     z                                                     *
C.    *  - E     energy of the particle at this point                  *
C.    *  - ELOS  the energy deposited into this block                  *
C.    * Example of one hit in that scheme:                             *
C.    *       EPHI 12                                                  *
C.    *       EZRI 41                                                  *
C.    *       BLOC  3                                                  *
C.    *        X    7.89 cm                                            *
C.    *        Y    -345.6 cm                                          *
C.    *        Z    1234.8 cm                                          *
C.    *        E    12 Gev                                             *
C.    *        ELOS  11.85 Gev                                         *
C.    * The FORTRAN  coding to  define the  set/det/hits information   *
C.    * could be:                                                      *
C.    *   DIMENSION NAMESV(3),NBITSV(3)                                *
C.    *   DIMENSION NAMESH(5),NBITSH(5),ORIG(5),FACT(5)                *
C.    *   DATA NAMESV/'EPHI','EZRI','BLOC'/                            *
C.    *   DATA NBITSV/6,6,3/                                           *
C.    *   DATA NAMESH/'X   ','Y   ','Z   ','E   ','ELOS'/              *
C.    *   DATA NBITSH/5*16/                                            *
C.    *   DATA ORIG/3*1000.,0.,0./                                     *
C.    *   DATA FACT/3*10.,2*100./                                      *
C.    *       CALL GSDET ('ECAL','BLOC',3,NAMESV,NBITSV,2,100,100,     *
C.    *      +                         ISET,IDET)                      *
C.    *       CALL GSDETH('ECAL','BLOC',5,NAMESH,NBITSH,ORIG,FACT)     *
C.    *   Returns  the hit  parameters  for  detector IUDET  of  set   *
C.    * IUSET.  All arguments are explained above.                     *
C..   *                                                                *
C.    *                                                                *
C.    *       JS = LQ(JSET-ISET)                                       *
C.    *       JD = LQ(JS-IDET)                                         *
C.    *       JDH= LQ(JD-1)                                            *
C.    *       IQ(JDH+4*I-3)= NAMESH(I)                                 *
C.    *       IQ(JDH+4*I-2)= NBITSH(I)                                 *
C.    *        Q(JDH+4*I-1)= ORIG(I)                                   *
C.    *        Q(JDH+4*I)  = FACT(I)                                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGEOM                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcmzfo.inc"
      PARAMETER (NHEMX=100)
      DIMENSION NBITSH(1),ORIG(1),FACT(1),NAMH(NHEMX)
      CHARACTER*4 NAMESH(1),IUSET,IUDET
      EQUIVALENCE (WS(1),NAMH(1))
C.
C.    ------------------------------------------------------------------
C.
      IF(JSET.LE.0)GO TO 90
      NSET=IQ(JSET-1)
      IF(NSET.LE.0)GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF(ISET.LE.0)GO TO 90
      JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.LE.0)GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF(IDET.LE.0)GO TO 90
      JD=LQ(JS-IDET)
C
      CALL MZBOOK(IXCONS,JDH,JD,-1,'SJDH',0,0,4*NH,IOSJDH,0)
C
      NW=0
      IF(NH.GT.0)THEN
         CALL UCTOH(NAMESH,NAMH,4,4*NH)
         K=32
         DO 30 I=1,NH
            NB=NBITSH(I)
            IF(NB.LT.0)NB=0
            IF(NB.GE.32)NB=0
            IQ(JDH+4*I-3)=NAMH(I)
            IQ(JDH+4*I-2)=NB
            Q(JDH+4*I-1)=ORIG(I)
            Q(JDH+4*I )=FACT(I)
            IF(FACT(I).LE.0.)Q(JDH+4*I)=1.
            IF(NB.LE.0)THEN
               NW=NW+1
               K=32
            ELSE
               K=K+NB
               IF(K.GT.32)THEN
                  K=NB
                  NW=NW+1
               ENDIF
            ENDIF
  30     CONTINUE
      ENDIF
      IQ(JD+3)=NW
      IQ(JD+4)=NH
      GO TO 99
C
C              Error
C
  90  WRITE(CHMAIL,1000)IUSET,IUDET
      CALL GMAIL(0,0)
 1000 FORMAT(' ***** GSDETH ERROR FOR SET ',A4,' OR DETECTOR ',A4)
C
  99  RETURN
      END

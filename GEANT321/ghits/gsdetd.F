*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GSDETD(IUSET,IUDET,ND,NAMESD,NBITSD)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Handling sensitive DETector Digitisation parameters        *
C.    *     ---------------------------------------------------        *
C.    *                                                                *
C.    * Defines digitisation  parameters for  detector IUDET  of set   *
C.    * IUSET.                                                         *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user detector identifier                             *
C.    * ND        number of elements per digitisation                  *
C.    * NAMESD    the  ND   variable  names  for   the  digitisation   *
C.    *           elements                                             *
C.    * NBITSD    the  ND  bit  numbers  for  packing  the  variable   *
C.    *           values.                                              *
C.    *   The  routine  is  used at  initialisation  time  once  the   *
C.    * geometrical  volumes  have  been  defined  to  describe  the   *
C.    * digitisation elements  and the way  to do packing  in memory   *
C.    * and on tape.  Let us use the same example as in GSDETH.  The   *
C.    * non  geometrical  information  we want  to  store  for  each   *
C.    * digitisation is for example:                                   *
C.    *  - ADC  pulse height in a lead glass block.                    *
C.    * Example of one digitisation in that scheme:                    *
C.    *       EPHI 12                                                  *
C.    *       EZRI 41                                                  *
C.    *       BLOC  3                                                  *
C.    *       ADC  789                                                 *
C.    * The FORTRAN  coding to  define the  digitisation information   *
C.    * could be:                                                      *
C.    *   DATA NAMESD/'ADC '/                                          *
C.    *   DATA NBITSD/16/                                              *
C.    *       CALL GSDETD('ECAL','BLOC',1,NAMESD,NBITSD)               *
C.    *   Returns the digitisation parameters  for detector IUDET of   *
C.    * set IUSET.  All arguments as explained in GSDETD.              *
C..   *                                                                *
C.    *       JS = LQ(JSET-ISET)                                       *
C.    *       JD = LQ(JS-IDET)                                         *
C.    *       JDD= LQ(JD-2)                                            *
C.    *       IQ(JDD+2*I-1)=NAMESD(I)                                  *
C.    *       IQ(JDD+2*I)  =NBITSD(I)                                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGEOM                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcmzfo.inc"
      PARAMETER (NDEMX=100)
      DIMENSION NBITSD(1),NAMD(NDEMX)
      CHARACTER*4 NAMESD(1),IUSET,IUDET
      EQUIVALENCE (WS(1),NAMD(1))
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
      CALL MZBOOK(IXCONS,JDD,JD,-2,'SJDD',0,0,2*ND,IOSJDD,0)
C
      NW=0
      IF(ND.GT.0)THEN
         CALL UCTOH(NAMESD,NAMD,4,4*ND)
         K=32
         DO 30 I=1,ND
            NB=NBITSD(I)
            IF(NB.LT.0)NB=0
            IF(NB.GE.32)NB=0
            IQ(JDD+2*I-1)=NAMD(I)
            IQ(JDD+2*I )=NB
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
      IQ(JD+5)=NW
      IQ(JD+6)=ND
      GO TO 99
C
C              Error
C
  90  WRITE(CHMAIL,1000)IUSET,IUDET
      CALL GMAIL(0,0)
 1000 FORMAT(' ***** GSDETD ERROR FOR SET ',A4,' OR DETECTOR ',A4)
C
  99  RETURN
      END

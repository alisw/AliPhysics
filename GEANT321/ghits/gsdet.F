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
      SUBROUTINE GSDET(IUSET,IUDET,NV,NAMESV,NBITSV,IDTYPE
     +                ,NWHI,NWDI,ISET,IDET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Defines volume parameters for detector IUDET of set IUSET*
C.    *                                                                *
C.    *             Input parameters                                   *
C.    * IUSET     set identifier (4 characters), user defined          *
C.    * IUDET     detector identifier  (4 characters),   name of  an   *
C.    *           existing volume                                      *
C.    * NV        number of volume descriptors                         *
C.    * NAMESV    vector of NV volume descriptors (4 characters)       *
C.    * NBITSV    vector of  NV bit numbers  for packing  the volume   *
C.    *           numbers                                              *
C.    * IDTYPE    detector type, user defined                          *
C.    * NWHI      number of words for the primary allocation of HITS   *
C.    *           banks                                                *
C.    * NWDI      number of words for the primary allocation of DIGI   *
C.    *           banks when first allocation not sufficient           *
C.    *                                                                *
C.    *             Output parameters                                  *
C.    * ISET      position of set in bank JSET                         *
C.    * IDET      position of detector in bank JS=IB(JSET-ISET)        *
C.    *              If ISET=0 or IDET=0  error                        *
C.    * Remarks:                                                       *
C.    * - The vector NAMESV (length NV)  contains the list of volume   *
C.    *   names which  permit the  identification of  every physical   *
C.    *   detector with detector name IUDET.    [See example in HITS   *
C.    *   110].                                                        *
C.    * - Each  element of  the vector  NBITSV (length  NV)  is  the   *
C.    *   number  of  bits  used  for  packing  the  number  of  the   *
C.    *   corresponding volume,  when building the packed identifier   *
C.    *   of a given physical detector.                                *
C.    * - For more details see the example given in GSDETH.            *
C.    * - The detector type  IDTYPE is not used  internally by GEANT   *
C.    *   and  can be  defined by  the user  to distinguish  quickly   *
C.    *   between various kinds of detectors,  in the routine GUSTEP   *
C.    *   for example.                                                 *
C.    *                                                                *
C.    *       IQ(JSET+ISET) = IUSET                                    *
C.    *       JS = LQ(JSET-ISET) + pointer to set parameters           *
C.    *       IQ(JS+IDET)=IUDET                                        *
C.    *       JD= LQ(JS-1)  = pointer to detector IDET                 *
C.    *       IQ(JD+1)=Total number of words to store packed volumes   *
C.    *       IQ(JD+2)=NV                                              *
C.    *       IQ(JD+3)=Number of words required per hit                *
C.    *       IQ(JD+4)=Number of different hits types                  *
C.    *       IQ(JD+5)=Number of words required per digit              *
C.    *       IQ(JD+6)=Number of different digit types                 *
C.    *       IQ(JD+7)=NWHI                                            *
C.    *       IQ(JD+8)=NWDI                                            *
C.    *       IQ(JD+9)=Number of paths through the JVOLUM tree         *
C.    *       IQ(JD+10)= For an alias only, IDET of main detector      *
C.    *       IQ(JD+2*I+9) = name of volume i = NAMESV(I)              *
C.    *       IQ(JD+2*I+10)= number of bits/volume = NBITSV(I)         *
C.    *                                                                *
C.    *            The Detector Set data structure JSET                *
C.    *            ------------------------------------                *
C.    *                                                                *
C.    *                                        | JSET                  *
C.    *    NSET            ISET                v         NSET          *
C.    *     ................................................           *
C.    *     |              | |               |  | Set names|           *
C.    *     ................................................           *
C.    *                     | JS                                       *
C.    *                     |                                          *
C.    *    NDET       IDET  v                    NDET                  *
C.    *     ........................................                   *
C.    *     |        |  |  | | Detector names      |                   *
C.    *     ........................................                   *
C.    *               | JD                                             *
C.    *      -3 -2 -1 v                                                *
C.    *     ................................................           *
C.    *     |  |  |  |  |   Parameters of GSDET            |           *
C.    *     ................................................           *
C.    *      |  |  |                                                   *
C.    *      |  |  |  JDH                                              *
C.    *      |  |  |                                                   *
C.    *      |  |  |           .............................           *
C.    *      |  |  ............| Parameters of GSDETH      |           *
C.    *      |  |              .............................           *
C.    *      |  |                                                      *
C.    *      |  | JDD                                                  *
C.    *      |  |                                                      *
C.    *      |  |              .............................           *
C.    *      |  ...............|  Parameters of GSDETD     |           *
C.    *      |                 .............................           *
C.    *      |                                                         *
C.    *      |  JDU                                                    *
C.    *      |                 .............................           *
C.    *      ..................| Parameters of GSDETU      |           *
C.    *                        .............................           *
C.    *  JS = LQ(JSET-ISET) pointer to detector set number ISET        *
C.    * The JSET data structure is filled by GSDET, GSDETH,  GSDETD,   *
C.    * GSDETU and eventually by GSDETA.                               *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGEOM                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcmzfo.inc"
      PARAMETER (NVMAX=20)
      DIMENSION NBITSV(1),NAV(NVMAX)
      CHARACTER*4 NAMESV,IUSET,IUDET
      EQUIVALENCE (WS(1),NAV(1))
C.
C.    ------------------------------------------------------------------
C.
      ISET=0
      IDET=0
      IF(NV.GT.15)GO TO 94
C
C             Check if volume IUDET has been defined
C
      IF(JVOLUM.LE.0)GO TO 90
      NVOLUM=IQ(JVOLUM-1)
      CALL GLOOK(IUDET,IQ(JVOLUM+1),NVOLUM,IVOL)
      IF(IVOL.EQ.0)GO TO 90
C
      IF(JSET.EQ.0)THEN
C
C               Create mother JSET bank
C
         CALL MZBOOK(IXCONS,JSET,JSET,1,'SETS',0,0,0,5,0)
         IQ(JSET-5)=0
         NSET=0
      ELSE
         NSET=IQ(JSET-1)
         CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
         IF (ISET.NE.0) GO TO 30
      ENDIF
C
C               Create JS bank
C
      CALL MZPUSH(IXCONS,JSET,1,1,'I')
      NSET=NSET+1
C
      ISET=NSET
      CALL UCTOH(IUSET,IQ(JSET+ISET),4,4)
      CALL MZBOOK(IXCONS,JS,JSET,-ISET,'SETS',0,0,0,5,0)
C
C            Check if detector has already been defined
C
  30  JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.NE.0)THEN
         CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
         IF(IDET.NE.0) GO TO 92
      ENDIF
C
C                If not, create detector bank
C
      CALL MZPUSH(IXCONS,JS,1,1,'I')
      NDET=NDET+1
      IDET=NDET
      CALL UCTOH(IUDET,IQ(JS+IDET),4,4)
      CALL MZBOOK(IXCONS,JD,JS,-IDET,'SEJD',4,4,10+2*NV,IOSEJD,0)
      IQ(JD-5)=10*ISET+IDET
C
      NW=0
      IF(NV.GT.0)THEN
         CALL UCTOH(NAMESV,NAV,4,4*NV)
         K=32
         DO 70 I=1,NV
            NB=NBITSV(I)
            IF(NB.LT.0)NB=0
            IF(NB.GE.32)NB=0
            IQ(JD+2*I+ 9)=NAV(I)
            IQ(JD+2*I+10)=NB
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
  70     CONTINUE
      ENDIF
  80  IQ(JD+1)=NW
      IQ(JD+2)=NV
      IQ(JD+7)=NWHI
      IQ(JD+8)=NWDI
      IQ(JD+10)=0
C
C             Now enter set,det into JVOLUM data structure
C
      CALL GSATT(IUDET,'SET ',ISET)
      CALL GSATT(IUDET,'DET ',IDET)
      CALL GSATT(IUDET,'DTYP',IDTYPE)
      GO TO 99
C
C              Errors
C
  90  WRITE(CHMAIL,1000)IUDET
      CALL GMAIL(0,0)
 1000 FORMAT(' ***** GSDET ERROR, VOLUME ',A4,' NOT DEFINED')
      GO TO 99
C
  92  WRITE(CHMAIL,2000)IUSET,IUDET
      CALL GMAIL(0,0)
 2000 FORMAT(' ***** GSDET ERROR ,SET ',A4, ' DETECTOR ',A4,
     +        ' ALREADY DEFINED')
      GO TO 99
C
  94  WRITE(CHMAIL,3000)NV
      CALL GMAIL(0,0)
 3000 FORMAT(' ***** GSDET ERROR ,SET ',A4, ' DETECTOR ',A4,
     +        ' Too many descriptors:',I5)
C
  99  RETURN
      END

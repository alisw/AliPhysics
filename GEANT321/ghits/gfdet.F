*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFDET(IUSET,IUDET,NV,NAMESV,NBITSV,IDTYPE
     +                ,NWHI,NWDI,ISET,IDET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       returns volume parameters for detector IUDET of set IUSET*
C.    *                                                                *
C.    *             Input parameters                                   *
C.    * IUSET     set identifier (4 characters), user defined          *
C.    * IUDET     detector identifier  (4 characters),   name of  an   *
C.    *           existing volume                                      *
C.    *                                                                *
C.    *             Output parameters                                  *
C.    * NV        number of volume descriptors                         *
C.    * NAMESV    vector of NV volume descriptors (4 characters)       *
C.    * NBITSV    vector of  NV bit numbers  for packing  the volume   *
C.    *           numbers                                              *
C.    * IDTYPE    detector type, user defined                          *
C.    * NWHI      number of words for the primary allocation of HITS   *
C.    *           banks                                                *
C.    * NWDI      number of words for the primary allocation of DIGI   *
C.    *           banks when first allocation not sufficient           *
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
C.    * - Vectors NAMESV and NBITSV must be dimensionned at least      *
C.    *   to NV in the calling routine.                                *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun , M.Maire *********                     *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION NBITSV(1)
      CHARACTER*4 NAMESV(1),IUSET,IUDET
C.
C.    ------------------------------------------------------------------
C.
      ISET=0
      IDET=0
C
C             Check if detector IUDET has been defined
C
      IF (JSET.LE.0) GO TO 90
      NSET = IQ(JSET-1)
      IF (NSET.LE.0) GO TO 90
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF (ISET.EQ.0) GO TO 90
      JS = LQ(JSET-ISET)
      NDET = IQ(JS-1)
      IF (NDET.LE.0) GO TO 90
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF (IDET.EQ.0) GO TO 95
      JD=LQ(JS-IDET)
      NV=IQ(JD+2)
      NWHI=IQ(JD+7)
      NWDI=IQ(JD+8)
C
      IF(NV.GT.0)THEN
         DO 10 I=1,NV
            CALL UHTOC(IQ(JD+2*I+ 9),4,NAMESV(I),4)
            NBITSV(I)=IQ(JD+2*I+10)
  10     CONTINUE
      ENDIF
C
      CALL GFATT(IUDET,'DTYP',IDTYPE)
      GO TO 99
C
   90 WRITE (CHMAIL, 1000) IUSET
      CALL GMAIL(0,0)
 1000 FORMAT (' ***** GFDET ERROR SET ',A4,' NOT FOUND')
      GO TO 99
   95 WRITE (CHMAIL, 2000) IUSET,IUDET
      CALL GMAIL(0,0)
 2000 FORMAT (' ***** GFDET ERROR FOR SET ',A4,
     + ' DETECTOR ',A4,' NOT FOUND')
C
  99  RETURN
      END

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
      SUBROUTINE GFDIGI(IUSET,IUDET,NTDIM,NVDIM,NDDIM,NDMAX,NUMVS
     +,                 LTRA,NTRA,NUMBV,KDIGI,NDIGS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *   Returns  the   digitisations  for   the  physical   volume   *
C.    * specified by the  list NUMVS with generic  volume name IUDET   *
C.    * belonging to set IUSET.                                        *
C.    * IUSET     user set identifier                                  *
C.    * IUDET     user    detector   identifier    (name   of    the   *
C.    *           corresponding sensitive volume)                      *
C.    * NTDIM     1st  dimension of  LTRA  (max.   number of  tracks   *
C.    *           contributing)                                        *
C.    * NVDIM     1st  dimension of  NUMVS,  NUMBV  (usually =NV, the  *
C.    *           number of volume descriptors which permit to iden-   *
C.    *           tify a given detector, possibly smaller than NV)     *
C.    * NDDIM     1st dimension of KDIGI (argument ND of GSDETD)       *
C.    * NDMAX     is  the  maximum  number of  digitisations  to  be   *
C.    *           returned                                             *
C.    * NUMVS     is a  1-Dim array that  must contain on  input the   *
C.    *           geometric  path  of  the  detector  volume  to  be   *
C.    *           selected.                                            *
C.    *           All 0  interpreted as  'all physical  volumes with   *
C.    *           generic name IUDET'                                  *
C.    * LTRA      is a 2-Dim  array that will contain  on output for   *
C.    *           each digitisation the numbers  of the tracks which   *
C.    *           have produced it                                     *
C.    * NTRA      is a 1-Dim  array that will contain  on output for   *
C.    *           each  digitisation  the  total  number  of  tracks   *
C.    *           contributing.                                        *
C.    *           In case this  number is greater than  NTDIM,  only   *
C.    *           the  first  NTDIM  corresponding   tracks  can  be   *
C.    *           returned on LTRA                                     *
C.    * NUMBV     is a 2-Dim  array that will contain  on output for   *
C.    *           each digitisation the list of volume numbers which   *
C.    *           identify each physical volume                        *
C.    * KDIGI     is  a 2-Dim  array  that  will contain  the  NDIGI   *
C.    *           digitisations                                        *
C.    * NDIGI     is  the  total  number of  digitisations  in  this   *
C.    *           detector.                                            *
C.    *           In  case  the  total number  of  digitisations  is   *
C.    *           greater than NDMAX,   NDIGI is set to  NDMAX+1 and   *
C.    *           only NDMAX digitisations are returned                *
C.    *     -  KDIGI(1,I)  =  digitisation type  1 for  digitisation   *
C.    *   number I                                                     *
C.    *     - NUMBV(1,I) = volume number 1 for digitisation number I   *
C.    *     -  LTRA (1,I)   =  first  track number  contributing  to   *
C.    *   digitisation number I                                        *
C.    *   In the calling  routine,  the arrays LTRA,   NTRA,  NUMVS,   *
C.    * NUMBV, KDIGI must be dimensioned to:                           *
C.    *    LTRA (NTDIM,NDMAX)                                          *
C.    *    NTRA (NDMAX)                                                *
C.    *    NUMVS(NVDIM)                                                *
C.    *    NUMBV(NVDIM,NDMAX)                                          *
C.    *    KDIGI(NDDIM,NDMAX)                                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    W.Gebel  *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      PARAMETER (NVMAX=20)
      DIMENSION NUMVT(NVMAX),NUMVS(NVDIM),NUMBV(NVDIM,1)
      DIMENSION LTRA(NTDIM,1),NTRA(1),KDIGI(NDDIM,1)
      EQUIVALENCE (WS(1),NUMVT(1))
      CHARACTER*4 IUSET,IUDET
C.
C.    --------------------------------------------------------------------------
C.
C             Find if selected set, detector exists
C
      NDIGS=0
      IF(JDIGI.LE.0)GO TO 999
      NSET=IQ(JSET-1)
      CALL GLOOK(IUSET,IQ(JSET+1),NSET,ISET)
      IF(ISET.LE.0)GO TO 999
C
      JS=LQ(JSET-ISET)
      JDI=LQ(JDIGI-ISET)
      IF(JS.LE.0)GO TO 999
      IF(JDI.LE.0)GO TO 999
      NDET=IQ(JS-1)
      CALL GLOOK(IUDET,IQ(JS+1),NDET,IDET)
      IF(IDET.EQ.0)GO TO 999
C
      JD=LQ(JS-IDET)
      JDID=LQ(JDI-IDET)
      IF(JDID.LE.0)GO TO 999
      JDDI=LQ(JD-2)
C
      ILAST=IQ(JDI+IDET)
      IF(ILAST.EQ.0)GO TO 999
      NV=IQ(JD+2)
      ND=IQ(JD+6)
C
C
C             Loop on all digits
C
C
      IDIG=0
      I=0
      NWDI=0
C
   10 CONTINUE
      I=I+NWDI
      IF(I.GE.ILAST)GO TO 110
      NWDI=IQ(JDID+I+1)
      NK=2
C
      NTRM1= IBITS(IQ(JDID+I+NK),0,16)
      NTRT = NTRM1+1
      NWTR = NTRT/2+1
      NK   = NK+NWTR
C
C             Find the selected volume
C             (if NO volumes exist take ALL digits)
C
      IF(NV.GT.0)THEN
         K=1
         DO 40 IV=1,NV
            NB=IQ(JD+2*IV+10)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   K=1
                   NK=NK+1
               ENDIF
               IF(IV.LE.NVMAX)NUMVT(IV)=IQ(JDID+I+NK)
               IF(IV.NE.NV)NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               IF(IV.LE.NVMAX)NUMVT(IV)=IBITS(IQ(JDID+I+NK),K-1,NB)
               K=K+NB
            ENDIF
            IF(IV.LE.NVDIM)THEN
               IF(NUMVS(IV).NE.0)THEN
                  IF(NUMVS(IV).NE.NUMVT(IV))GO TO 10
               ENDIF
            ENDIF
   40    CONTINUE
         NK=NK+1
      ENDIF
C
C
C
C ========>   Now store number of tracks and volume numbers,
C             and fetch track numbers and digits
C
      IDIG=IDIG+1
      IF(IDIG.GT.NDMAX)GO TO 110
C
      NTRA(IDIG)=NTRT
      NVMIN=MIN(NV,NVDIM)
      CALL VZERO (NUMBV(1,IDIG),NVDIM)
      CALL UCOPY (NUMVT(1),NUMBV(1,IDIG),NVMIN)
C
C             Get track numbers
C
      MK=NK
      NK=2
      IF(NTRT.GT.0)THEN
         IF(NTRM1.GE.1)THEN
            DO 54 ITR=1,NTRM1,2
               IF(ITR.LE.NTDIM)THEN
                  LTRA(ITR ,IDIG)=IBITS(IQ(JDID+I+NK),16,16)
               ENDIF
               NK=NK+1
               IF(ITR.LT.NTDIM)THEN
                  LTRA(ITR+1,IDIG)=IBITS(IQ(JDID+I+NK), 0,16)
               ENDIF
   54       CONTINUE
         ENDIF
         IF(NTRT.LE.NTDIM)THEN
           IF(MOD(NTRT,2).EQ.1)
     +               LTRA(NTRT,IDIG)=IBITS(IQ(JDID+I+NK),16,16)
         ENDIF
      ENDIF
      NK=MK
C
C             Get unpacked digits
C
      IF(ND.LE.0)GO TO 10
      K=1
      DO 90 ID=1,ND
         NB=IQ(JDDI+2*ID)
         IF(NB.LE.0)THEN
            IF(K.GT.1)THEN
                K=1
                NK=NK+1
            ENDIF
            IF(ID.LE.NDDIM)KDIGI(ID,IDIG)=IQ(JDID+I+NK)
            IF(ID.NE.ND)NK=NK+1
         ELSE
            IF(K+NB.GT.33)THEN
               K=1
               NK=NK+1
            ENDIF
            IF(ID.LE.NDDIM)KDIGI(ID,IDIG)=IBITS(IQ(JDID+I+NK),K-1,NB)
            K=K+NB
         ENDIF
  90  CONTINUE
C
      GO TO 10
C
 110  NDIGS=IDIG
C
 999  RETURN
      END

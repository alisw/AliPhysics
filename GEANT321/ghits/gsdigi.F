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
      SUBROUTINE GSDIGI(ISET,IDET,LTRA,NTRA,NUMBV,KDIGI,IDIG)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Routines to Communicate with the data structure JDIGI       *
C.    *   ---------------------------------------------------------    *
C.    *                                                                *
C.    *   Stores element  values for  current digitisation  into the   *
C.    * data structure JDIGI.                                          *
C.    * ISET      set number                                           *
C.    * IDET      detector number                                      *
C.    * LTRA      list  of   NTRA  track   numbers  producing   this   *
C.    *           digitisation                                         *
C.    * NUMBV     volume  numbers corresponding  to  list NAMESV  of   *
C.    *           GSDET                                                *
C.    * KDIGI     integer  array of  values  for current  digisation   *
C.    *           elements                                             *
C.    * IDIG      on return, current digitisation number.              *
C.    *           If =0 digitisation has not been stored.              *
C.    *                                                                *
C.    *           The Digitisation data structure JDIGI                *
C.    *           -------------------------------------                *
C.    *                                                                *
C.    *                                            | JDIGI             *
C.    *    NSET                  ISET              v                   *
C.    *     ..........................................                 *
C.    *     |                   |  |              |  |                 *
C.    *     ..........................................                 *
C.    *                           |                                    *
C.    *                           | JDIG                               *
C.    *      NDET    IDET         v               NDET                 *
C.    *        .....................................                   *
C.    *        |      |  |      |  |               |                   *
C.    *        .....................................                   *
C.    *                |                                               *
C.    *                | JDD                                           *
C.    *                v                                               *
C.    *               ...........................................      *
C.    *               |  | 1st digitisation | 2nd digitisation, etc.   *
C.    *               ...........................................      *
C.    *                         Bank layout                            *
C.    * JDIG           = LQ(JDIGI-ISET),   pointer to  digitisations   *
C.    *                for set ISET                                    *
C.    * JDD            = LQ(JDIG-IDET),  pointer to digitisations of   *
C.    *                detector IDET of set ISET                       *
C.    * IQ(JDIG+IDET)  pointer to last word of last digitisation for   *
C.    *                this detector                                   *
C.    * IQ(JDD+1)      1st word of first digitisation                  *
C.    * IQ(JDD+NWD+1)  1st word of second digitisation                 *
C.    *                JS=LQ(JSET-ISET)                                *
C.    *                JD=LQ(JS-IDET)                                  *
C.    *                NWD=IQ(JD+5)                                    *
C.    *   The JDIGI  structure is  filled with  the routine  GSDIGI.   *
C.    * The routine GFDIGI can be used  to get the digitisations for   *
C.    * a detector IDET and set ISET.                                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *       Authors    R.Brun, W.Gebel  *********                    *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      COMMON/GCLOCA/NLOCAL(2),JS,JD,JDDI,JDI,JDID,LOCAL(15)
      DIMENSION NUMBV(1),KDIGI(1),LTRA(1)
C.
C.    ------------------------------------------------------------------
C.
C             Find if selected set, detector exists
C
      IDIG=0
C
      IF(JSET.LE.0)GO TO 99
      NSET=IQ(JSET-1)
      IF(NSET.LE.0)GO TO 99
      IF(ISET.LE.0)GO TO 99
      IF(ISET.GT.NSET)GO TO 99
C
      JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.LE.0)GO TO 99
      IF(IDET.LE.0)GO TO 99
      IF(IDET.GT.NDET)GO TO 99
C
      JD=LQ(JS-IDET)
      JDDI=LQ(JD-2)
      IF(JDDI.LE.0)GO TO 99
      NW=IQ(JD+1)+IQ(JD+5)+2
      ND=IQ(JD+6)
      NV=IQ(JD+2)
      NWTR=NTRA/2+1
      NWD=NW+NWTR
C
C              Create DIGItisation master bank
C
      IF(JDIGI.EQ.0)THEN
         CALL MZBOOK(IXDIV,JDIGI,JDIGI,1,'DIGI',NSET,NSET,0,2,0)
         IQ(JDIGI-5)=0
      ENDIF
      JDI=LQ(JDIGI-ISET)
      IF(JDI.EQ.0)THEN
         CALL MZBOOK(IXDIV,JDI,JDIGI,-ISET,'DIGI',NDET,NDET,NDET,2,0)
      ENDIF
C
      JDID=LQ(JDI-IDET)
      IF(JDID.EQ.0)THEN
C
C            Create DIGItisation bank
C
         NWDI=IQ(JD+8)
         CALL MZBOOK(IXDIV,JDID,JDI,-IDET,'SJDI',0,0,NWDI,1,0)
         IQ(JDI+IDET)=0
      ENDIF
C
C           Check if enough space. If not increase bank size
C
  10  NDID=IQ(JDID-1)
      ILAST=IQ(JDI+IDET)
      NFREE=NDID-ILAST
      IF(NFREE.LE.NWD)THEN
         NWDI2=MAX(100,NWD,IQ(JD+8)/2)
         CALL MZPUSH(IXDIV,JDID,0,NWDI2,'I')
         GO TO 10
      ENDIF
C
      IQ(JDI+IDET)=ILAST+NWD
      IF(ILAST.NE.0)IDIG=IQ(JDID+ILAST)
      IDIG=IDIG+1
C
C
C
C
C ========>    Store tracks numbers, volumes numbers and digits
C
      IQ(JDID+ILAST+1)=NWD
      IQ(JDID+ILAST+NWD)=IDIG
      NK=ILAST+2
C
C             Store packed track numbers
C             Every 2 consecutive numbers into 1 word
C             1st half of 1st word: NTRA-1
C
      IQ(JDID+NK)=0
      IF(NTRA.GT.0)THEN
         NTRM1=NTRA-1
         CALL MVBITS(NTRM1,0,16,IQ(JDID+NK),0)
         IF(NTRM1.GE.1)THEN
            DO 23 ITR=1,NTRM1,2
               CALL MVBITS(LTRA(ITR),0,16,IQ(JDID+NK),16)
               NK=NK+1
               CALL MVBITS(LTRA(ITR+1),0,16,IQ(JDID+NK),0)
   23       CONTINUE
         ENDIF
         IF(MOD(NTRA,2).EQ.1)
     +          CALL MVBITS(LTRA(NTRA),0,16,IQ(JDID+NK),16)
      ENDIF
      NK=NK+1
C
C           Store packed volume numbers
C
      IF(NV.GT.0)THEN
         K=1
C
         DO 50 I=1,NV
            NB=IQ(JD+2*I+10)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               IQ(JDID+NK)=NUMBV(I)
                K=1
               IF(I.NE.NV)NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               CALL MVBITS(NUMBV(I),0,NB,IQ(JDID+NK),K-1)
               K=K+NB
            ENDIF
  50     CONTINUE
         NK=NK+1
      ENDIF
C
C          Store packed digits
C
      IF(ND.GT.0)THEN
         K=1
         DO 90 I=1,ND
            NB=IQ(JDDI+2*I)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               IQ(JDID+NK)=KDIGI(I)
               NK=NK+1
               K=1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               CALL MVBITS(KDIGI(I),0,NB,IQ(JDID+NK),K-1)
               K=K+NB
            ENDIF
  90     CONTINUE
      ENDIF
C
  99  RETURN
      END

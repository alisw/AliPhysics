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
      SUBROUTINE GSAHIT(ISET,IDET,ITRA,NUMBV,HITS,IHIT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Routines to Communicate with the JHITS data structure       *
C.    *   ---------------------------------------------------------    *
C.    *                                                                *
C.    *   Stores  element  values  for current  hit  into  the  data   *
C.    * structure JHITS.                                               *
C.    * ISET      set number  (can be obtained from  /GCSETS/ filled   *
C.    *           by GFINDS)                                           *
C.    * IDET      detector    number   "           "               "   *
C.    *           "                                                    *
C.    * ITRA      track number producing this hit                      *
C.    * NUMBV     array  of  volume numbers  corresponding  to  list   *
C.    *           NAMESV of GSDET                                      *
C.    * HITS      array of values for current hit elements             *
C.    * IHIT      on return, current hit number.  If =0, hit has not   *
C.    *           been stored.                                         *
C.    *                                                                *
C.    *       JH=LQ(JHITS-ISET)                                        *
C.    *       JHD=LQ(JH-IDET)                                          *
C.    *       IQ(JH+IDET)= pointer to LAST USED word in JHD            *
C.    *         Each hit is packed into JHD in the following format    *
C.    *        --Track number ITRA not packed                          *
C.    *        --Volume numbers packed                                 *
C.    *        --Hits transformed and packed                           *
C.    *                                                                *
C.    *                The Hit data structure JHITS                    *
C.    *                ----------------------------                    *
C.    *                                                                *
C.    *                                            | JHITS             *
C.    *    NSET                  ISET              v                   *
C.    *     ..........................................                 *
C.    *     |                   |  |              |  |                 *
C.    *     ..........................................                 *
C.    *                           |                                    *
C.    *                           | JH                                 *
C.    *      NDET    IDET         v               NDET                 *
C.    *        .....................................                   *
C.    *        |      |  |      |  |               |                   *
C.    *        .....................................                   *
C.    *                |                                               *
C.    *                | JHD                                           *
C.    *                v                                               *
C.    *               .........................................        *
C.    *               |  | 1st hit        | 2nd hit, etc.     |        *
C.    *               .........................................        *
C.    *                         Bank layout                            *
C.    * JH            =  LQ(JHITS-ISET,)  pointer  to  hits for  set   *
C.    *           number ISET                                          *
C.    * JHD           = LQ(JH-IDET),   pointer to  hits of  detector   *
C.    *           IDET                                                 *
C.    *                  of set ISET                                   *
C.    * IQ(JH+IDET)    number of words  used so far for  storing the   *
C.    *           hits                                                 *
C.    *                  of detector IDET                              *
C.    * IQ(JHD+1)       1st word of 1st hit                            *
C.    * IQ(JHD+NWH+1)  1st word of 2nd hit                             *
C.    *                      JS=LQ(JSET-ISET)                          *
C.    *                      JD=LQ(JS-IDET)                            *
C.    *                      NWH=IQ(JD+3)                              *
C.    *   The JHITS structure is filled with the routines GSAHIT and   *
C.    * GSCHIT.   The routine GFHITS can be used to get the hits for   *
C.    * a detector IDET and set ISET.                                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSTEP                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      COMMON/GCLOCA/NLOCAL(2),JS,JD,JDH,JH,JHD,LOCAL(15)
      DIMENSION NUMBV(1),HITS(1)
      SAVE NMESS
      DATA NMESS/0/
C.
C.    ------------------------------------------------------------------
C.
      IHIT=0
      IF(JSET.LE.0)GO TO 99
      NSET=IQ(JSET-1)
      IF(NSET.LE.0)GO TO 99
      IF(ISET.LE.0)GO TO 99
      IF(ISET.GT.NSET)GO TO 99
      JS=LQ(JSET-ISET)
      NDET=IQ(JS-1)
      IF(NDET.LE.0)GO TO 99
      IF(IDET.LE.0)GO TO 99
      IF(IDET.GT.NDET)GO TO 99
      JD=LQ(JS-IDET)
      JDH=LQ(JD-1)
      IF(JDH.LE.0)GO TO 99
      NW=IQ(JD+1)+IQ(JD+3)+1
      NH=IQ(JD+4)
      NV=IQ(JD+2)
C
C              Create HITS master bank
C
      IF(JHITS.EQ.0)THEN
         CALL MZBOOK(IXDIV,JHITS,JHITS,1,'HITS',NSET,NSET,0,2,0)
         IQ(JHITS-5)=0
      ENDIF
      JH=LQ(JHITS-ISET)
      IF(JH.EQ.0)THEN
         CALL MZBOOK(IXDIV,JH,JHITS,-ISET,'HITS',NDET,NDET,NDET,2,0)
      ENDIF
C
      JHD=LQ(JH-IDET)
      IF(JHD.EQ.0)THEN
C
C            Create Hits bank
C
         NWHI=IQ(JD+7)
         CALL MZBOOK(IXDIV,JHD,JH,-IDET,'SJHD',0,0,NWHI,1,0)
         IQ(JHD-5)=1000*ISET+IDET
         ILAST=0
      ELSE
C
C           Check if enough space. If not increase bank size
C
         NHD=IQ(JHD-1)
         ILAST=IQ(JH+IDET)
         NFREE=NHD-ILAST
         IF(NFREE.LE.NW)THEN
            NWHI2=MAX(100,NW,IQ(JD+7)/2)
            CALL MZPUSH(IXDIV,JHD,0,NWHI2,'I')
            JS  = LQ(JSET-ISET)
            JD  = LQ(JS-IDET)
            JDH = LQ(JD-1)
            JH  = LQ(JHITS-ISET)
         ENDIF
      ENDIF
C
      IQ(JH+IDET)=ILAST+NW
      IHIT=IQ(JH+IDET)/NW
C
C ========>    Store track number,volumes numbers and hits
C
      IQ(JHD+ILAST+1)=ITRA
C
C
C           Store packed volume numbers
C
      NK=ILAST+2
      IF(NV.GT.0)THEN
         K=1
C
         DO 50 I=1,NV
            NB=IQ(JD+2*I+10)
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               IQ(JHD+NK)=NUMBV(I)
               K=1
               IF(I.NE.NV)NK=NK+1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               CALL MVBITS(NUMBV(I),0,NB,IQ(JHD+NK),K-1)
               K=K+NB
            ENDIF
  50     CONTINUE
         NK=NK+1
      ENDIF
C
C          Store packed hits
C          Before packing, hits are changed to integers
C          Origin is shifted to have only positive integers
C          Result is multiplied by a constant to get resolution
C
      IF(NH.GT.0)THEN
         K=1
         DO 90 I=1,NH
            NB=IQ(JDH+4*I-2)
            XHIT=(HITS(I)+Q(JDH+4*I-1))*Q(JDH+4*I)
            IF(NB.EQ.0)THEN
               VALMX=2.147483E+9
            ELSE
               VALMX=2.**NB-1.
            ENDIF
            IFLAG=0
            IF(XHIT.LT.0.)THEN
               XHIT=0.
               IFLAG=1
            ELSE IF(XHIT.GT.VALMX)THEN
               XHIT=VALMX
               IFLAG=1
            ENDIF
            IF(IFLAG.NE.0)THEN
               NMESS=NMESS+1
               IF(NMESS.LT.10)THEN
                  WRITE(CHMAIL,1000)IQ(JSET+ISET),IQ(JS+IDET),I,HITS(I)
     +             ,Q(JDH+4*I-1),Q(JDH+4*I)
                  CALL GMAIL(0,0)
               ENDIF
            ENDIF
C
            KHIT=XHIT+0.5
            IF(NB.LE.0)THEN
               IF(K.GT.1)THEN
                   NK=NK+1
               ENDIF
               IQ(JHD+NK)=KHIT
               NK=NK+1
               K=1
            ELSE
               IF(K+NB.GT.33)THEN
                  K=1
                  NK=NK+1
               ENDIF
               CALL MVBITS(KHIT,0,NB,IQ(JHD+NK),K-1)
               K=K+NB
            ENDIF
  90     CONTINUE
      ENDIF
C
 1000 FORMAT(' ***** GSAHIT OVERFLOW WHEN IUSET= ',A4,' IUDET= ',
     +A4,' HITS(',I2,')=',E14.7,' ORIG= ',E14.7,' FACT= ',E14.7)
  99  RETURN
      END

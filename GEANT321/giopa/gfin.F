*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFIN(LUN,CHOBJ,NKEYS,IDVERS,CHOPT,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to read GEANT object(s) fromin the FZ file       *
C.    *       The data structures from disk are read in memory         *
C.    *                                                                *
C.    *       LUN    Logical unit                                      *
C.    *                                                                *
C.    *       CHOBJ  The type of data structure to be read:            *
C.    *              MATE material                                     *
C.    *              TMED tracking medium                              *
C.    *              VOLU volumes                                      *
C.    *              ROTM rotation matrix                              *
C.    *              SETS detector set                                 *
C.    *              PART particle                                     *
C.    *              SCAN geometry                                     *
C.    *              INIT read all above                               *
C.    *              KINE this keyword will trigger the read of        *
C.    *                   KINE and VERT unless the flag 'S' is set     *
C.    *              DIGI digitisation                                 *
C.    *              DRAW drawing                                      *
C.    *              HEAD event header                                 *
C.    *              HITS hits                                         *
C.    *              RUNG run                                          *
C.    *              STAK particle temporary stack                     *
C.    *              STAT volume statistic                             *
C.    *              VERT vertex                                       *
C.    *              JXYZ track points                                 *
C.    *              TRIG this keyword will trigger the read of        *
C.    *                   DIGI, HEAD, HITS, KINE, VERT abd JXYZ        *
C.    *                   unless the 'S' flag is set                   *
C.    *                                                                *
C.    *       NKEYS  number of keys in vector CHOBJ                    *
C.    *                                                                *
C.    *       IDVERS version of the data structure to be read in       *
C.    *                                                                *
C.    *       CHOPT  List of options                                   *
C.    *                   'I'      read only initialisation data       *
C.    *                            structures                          *
C.    *                   'K'      read only KINE and TRIG data        *
C.    *                            structures                          *
C.    *                   'T'      read only DIGI, HEAD, HITS, KINE,   *
C.    *                            VERT and JXYZ data structures       *
C.    *              even if other keys are specified in CHOBJ         *
C.    *                                                                *
C.    *                   'S'       interpret KINE to mean only        *
C.    *                             KINE and TRIG and INIT to mean     *
C.    *                             nothing                            *
C.    *                   'Q'       quiet option, no message is        *
C.    *                             printed                            *
C.    *                                                                *
C.    *       IER    error flag. <0 ZEBRA error flag as returned in    *
C.    *                             IQUEST(1)                          *
C.    *                           0 read completed successfully        *
C.    *                          >0 if only IER structures read in     *
C.    *                                                                *
C.    *    The FZ data base has been created via GOPEN/GFOUT           *
C.    *                                                                *
C.    *                                                                *
C.    *      Example.                                                  *
C.    *                                                                *
C.    *      CALL GOPEN(1,'I',1024,IER)                                *
C.    *      CALL GFIN (1,'VOLU',1,0,' ',IER)                          *
C.    *      CALL GFIN (1,'MATE',1,0,' ',IER)                          *
C.    *      CALL GFIN (1,'TMED',1,0,' ',IER)                          *
C.    *      CALL GFIN (1,'ROTM',1,0,' ',IER)                          *
C.    *      CALL GFIN (1,'PART',1,0,' ',IER)                          *
C.    *      CALL GFIN (1,'SCAN',1,0,' ',IER)                          *
C.    *      CALL GFIN (1,'SETS',1,0,' ',IER)                          *
C.    *                                                                *
C.    *    ==>Called by : <USER> ,GOPEN                                *
C.    *       Author    F.Carminati *******                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcnum.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctime.inc"
*      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
*     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
*     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
      COMMON/QUEST/IQUEST(100)
      PARAMETER (NLINIT=9,NLKINE=2,NLTRIG=6,NMKEY=22)
      DIMENSION JNAMES(20),LINIT(NLINIT),LKINE(NLKINE)
      DIMENSION LTRIG(NLTRIG),IXD(NMKEY)
      DIMENSION LINK(NMKEY),IVERSI(NMKEY),LDIV(2),IRESUL(NMKEY)
      DIMENSION IDOLD(8), IDNEW(8), VEROLD(8), VERNEW(8)
      DIMENSION IUHEAD(2),ITRAN(23)
      EQUIVALENCE (JNAMES(1),JDIGI)
      CHARACTER*4 KNAMES(NMKEY),CHOBJ(*)
      CHARACTER*(*) CHOPT
      DATA KNAMES/'DIGI','DRAW','HEAD','HITS','KINE','MATE','PART',
     +     'ROTM','RUNG','SETS','STAK','STAT','TMED','NULL','VERT',
     +     'VOLU','JXYZ','NULL','NULL','NULL','SCAN','NULL'/
      DATA ITRAN/7,6,13,16,8,10,2,9,8*0,3,15,5,17,4,1,21/
      DATA IXD/2,1,2,2,2,8*1,2,2,1,2,3*0,1,0/
      DATA LINIT/2,6,7,8,9,10,13,16,21/
      DATA LKINE/5,15/
      DATA LTRIG/1,3,4,5,15,17/
      DATA IDNEW / 8*0 /
      DATA VERNEW / 8*0. /
C.
C.    ------------------------------------------------------------------
C.
      IQUEST(1)=0
      LDIV(1)  =IXCONS
      LDIV(2)  =IXDIV
      KVOL=JVOLUM
      IER=0
*
      IOPTI=INDEX(CHOPT,'i')+INDEX(CHOPT,'I')
      IOPTT=INDEX(CHOPT,'t')+INDEX(CHOPT,'T')
      IOPTK=INDEX(CHOPT,'k')+INDEX(CHOPT,'K')
      IOPTS=INDEX(CHOPT,'s')+INDEX(CHOPT,'S')
      IOPTQ=INDEX(CHOPT,'q')+INDEX(CHOPT,'Q')
*
*     Save old JRUNG dates and versions
      IF(JRUNG.GT.0) THEN
         DO 10 J=1,8
            IDOLD(J) = IQ(JRUNG+10+J)
            VEROLD(J) = Q(JRUNG+20+J)
   10    CONTINUE
      ENDIF
*
      NLINK=0
      DO 100 JKEY=1,NKEYS
         IF(IOPTS.EQ.0) THEN
         IF(CHOBJ(JKEY).EQ.'INIT') THEN
            DO 30 J=1, NLINIT
               DO 20  MLINK=1,NLINK
                  IF(LINK(MLINK).EQ.LINIT(J)) GO TO 30
   20          CONTINUE
               NLINK=NLINK+1
               LINK(NLINK)=LINIT(J)
   30       CONTINUE
            GO TO 100
         ELSEIF(CHOBJ(JKEY).EQ.'TRIG') THEN
            DO 50 J=1, NLTRIG
               DO 40  MLINK=1,NLINK
                  IF(LINK(MLINK).EQ.LTRIG(J)) GO TO 50
   40          CONTINUE
               NLINK=NLINK+1
               LINK(NLINK)=LTRIG(J)
   50       CONTINUE
            GO TO 100
         ELSEIF(CHOBJ(JKEY).EQ.'KINE') THEN
            DO 70 J=1, NLKINE
               DO 60  MLINK=1,NLINK
                  IF(LINK(MLINK).EQ.LKINE(J)) GO TO 70
   60          CONTINUE
               NLINK=NLINK+1
               LINK(NLINK)=LKINE(J)
   70       CONTINUE
            GO TO 100
         ENDIF
         ENDIF
            DO 90 J=1,NMKEY
               IF(CHOBJ(JKEY).EQ.KNAMES(J)) THEN
                  DO 80 MLINK=1,NLINK
                     IF(LINK(MLINK).EQ.J) GO TO 100
   80             CONTINUE
                  NLINK=NLINK+1
                  LINK(NLINK)=J
                  GO TO 100
               ENDIF
   90       CONTINUE
            WRITE(CHMAIL,10300) CHOBJ(JKEY)
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
  100 CONTINUE
*
      IF(IOPTI.GT.0) THEN
         DO 120 J=1, NLINK
            DO 110 K=1, NLINIT
               IF(LINK(J).EQ.LINIT(K)) GO TO 120
  110       CONTINUE
            WRITE(CHMAIL,10000) KNAMES(LINK(J))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            LINK(J)=0
  120    CONTINUE
      ELSEIF(IOPTK.GT.0) THEN
         DO 140 J=1, NLINK
            DO 130 K=1, NLKINE
               IF(LINK(J).EQ.LKINE(K)) GO TO 140
  130       CONTINUE
            WRITE(CHMAIL,10100) KNAMES(LINK(J))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            LINK(J)=0
  140    CONTINUE
      ELSEIF(IOPTT.GT.0) THEN
         DO 160 J=1, NLINK
            DO 150 K=1, NLTRIG
               IF(LINK(J).EQ.LTRIG(K)) GO TO 160
  150       CONTINUE
            WRITE(CHMAIL,10200) KNAMES(LINK(J))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            LINK(J)=0
  160    CONTINUE
      ENDIF
      IOFF=0
      DO 170 J=1, NLINK
         IF(LINK(J).EQ.0) THEN
            IOFF=IOFF-1
         ELSE
            LINK(J+IOFF)=LINK(J)
         ENDIF
  170 CONTINUE
      NLINK=NLINK+IOFF
      NPOS=0
      DO 171 JL=1,NLINK
         IF(LINK(JL).EQ.9.OR.LINK(JL).EQ.3) THEN
            NPOS=JL
            GOTO 172
         ENDIF
  171 CONTINUE
  172 CONTINUE
*
      IF(NLINK.LE.0) THEN
         WRITE(CHMAIL,10400)
         IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         IER=-1
         GOTO 999
      ENDIF
*
      IF(IOPTI+IOPTK+IOPTT.EQ.0) THEN
*
*        We have to choose which event header to read, JRUNG or JHEAD
*        If the banks list contains banks depending from both headers,
*        the result is unpredictable. Error message to be inserted later.
         DO 168 J=1, NLINK
            DO 161 L=1, NLINIT
               IF(LINK(J).EQ.LINIT(L)) THEN
                  IOPTI=-1
                  GOTO 169
               ENDIF
  161       CONTINUE
            DO 162 L=1, NLKINE
               IF(LINK(J).EQ.LKINE(L)) THEN
                  IOPTK=-1
                  GOTO 169
               ENDIF
  162       CONTINUE
            DO 163 L=1, NLTRIG
               IF(LINK(J).EQ.LTRIG(L)) THEN
                  IOPTT=-1
                  GOTO 169
               ENDIF
  163       CONTINUE
  168    CONTINUE
  169    CONTINUE
      ENDIF
*
      DO 180 J=1, NLINK
         IVERSI(J)=0
         IRESUL(J)=0
  180 CONTINUE
*
*               Go for next start of event data structure
  190 IF(IOPTI.NE.0) THEN
         IF(JRUNG.NE.0)CALL MZDROP(IXCONS,JRUNG,' ')
         NUH=2
         CALL FZIN(LUN,IXCONS,JRUNG,1,'E',NUH,IUHEAD)
         IF(IQUEST(1).GE.2) THEN
            IER = -IQUEST(1)
            GO TO 240
         ENDIF
         IF(NPOS.NE.0) THEN
            IVERSI(NPOS)=IUHEAD(1)
            IRESUL(NPOS)=1
         ENDIF
      ELSEIF(IOPTT+IOPTK.NE.0) THEN
         IF(JHEAD.NE.0)CALL MZDROP(IXDIV,JHEAD,' ')
         NUH=2
         CALL FZIN(LUN,IXDIV,JHEAD,1,'E',NUH,IUHEAD)
         IF(IQUEST(1).GE.2) THEN
            IER = -IQUEST(1)
            GO TO 240
         ENDIF
         IF(NPOS.NE.0) THEN
            IVERSI(NPOS)=IUHEAD(1)
            IRESUL(NPOS)=1
         ENDIF
      ENDIF
*
      IVERIN = IUHEAD(1)
      IF(IDVERS.NE.0.AND.IDVERS.NE.IVERIN) THEN
         DO 200 I=1,NLINK
            LINK(I)=-ABS(LINK(I))
  200    CONTINUE
         GOTO 190
      ELSE
         IF (IDVERS .EQ. 0) IDVERS = IVERIN
         DO 210 I=1,NLINK
            LINK(I)= ABS(LINK(I))
  210    CONTINUE
      ENDIF
      NK   = IUHEAD(2)
      IF(NK.GT.10) THEN
         WRITE(CHMAIL,11100) NK
         IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         IER=-1
         GO TO 999
      ENDIF
      DO 230 IK=1,NK
C
C              Read next header
C
         NUH=2
         CALL FZIN(LUN,0,0,0,'S',NUH,IUHEAD)
         IF(IQUEST(1).GT.2)GO TO 320
         IKEY=ITRAN(IUHEAD(1))
         DO 220 I=1,NLINK
            NKEY=LINK(I)
            IF(IKEY.EQ.NKEY)THEN
               IDIV=LDIV(IXD(NKEY))
               IF(NKEY.LE.20)THEN
                  IF(JNAMES(NKEY).NE.0)THEN
                     CALL MZDROP(IDIV,JNAMES(NKEY),'L')
                     JNAMES(NKEY)=0
                  ENDIF
                  CALL FZIN(LUN,IDIV,JNAMES(NKEY),1,'A',NUH,IUHEAD)
               ELSE
                  NKL=NKEY-20
                  IF(ISLINK(NKL).NE.0)THEN
                     CALL MZDROP(IDIV,ISLINK(NKL),'L')
                     ISLINK(NKL)=0
                  ENDIF
                  CALL FZIN(LUN,IDIV,ISLINK(NKL),1,'A',NUH,IUHEAD)
               ENDIF
               IF(IQUEST(1).LE.2.AND.IQUEST(1).GE.0) THEN
                  IVERSI(I)=IVERIN
                  IRESUL(I)=1
                  GOTO 230
               ELSE
                  GOTO 320
               ENDIF
            ENDIF
  220    CONTINUE
  230 CONTINUE
*
  240 NIN=0
      DO 250 I=1,NLINK
         IF(IRESUL(I).EQ.1) THEN
            WRITE(CHMAIL,10500) KNAMES(LINK(I)), IVERSI(I)
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            NIN=NIN+1
         ELSEIF(LINK(I).GT.0) THEN
            WRITE(CHMAIL,10600) KNAMES(LINK(I))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         ELSEIF(LINK(I).LT.0) THEN
            WRITE(CHMAIL,10700) KNAMES(-LINK(I)), IDVERS
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         ENDIF
  250 CONTINUE
*
      IF(NIN.EQ.0) THEN
         WRITE(CHMAIL,10800)
         IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         IF(IER.GE.0) IER=-1
         GOTO 999
      ELSEIF(NIN.LT.NLINK) THEN
         IER=NIN
      ENDIF
*
      IF(KVOL.NE.JVOLUM)THEN
         NVOLUM=IQ(JVOLUM-1)
         CALL MZGARB(IXCONS,0)
         CALL GGDVLP
         CALL GGNLEV
      ENDIF
*
      IF(JVOLUM.GT.0) THEN
         NLEVEL=0
         NVOLUM=0
         DO 260 J=1, IQ(JVOLUM-2)
            IF(LQ(JVOLUM-J).EQ.0) GOTO 270
            NVOLUM=NVOLUM+1
  260    CONTINUE
  270    CONTINUE
      ENDIF
*
      IF(JTMED.NE.0 )THEN
         CALL UCOPY(Q(JTMED+1),CUTGAM,10)
         NTMED=IQ(JTMED-2)
      ENDIF
*
      IF(JPART.NE.0 ) NPART = IQ(JPART-2)
      IF(JVERTX.NE.0) NVERTX = IQ(JVERTX+1)
      IF(JKINE.NE.0) NTRACK = IQ(JKINE+1)
      IF(JMATE.NE.0 ) NMATE = IQ(JMATE-2)
      IF(JROTM.NE.0 ) NROTM = IQ(JROTM-2)
      IF(JDRAW.GT.0 ) THEN
         NKVIEW = IQ(JDRAW-2)
      ELSE
         NKVIEW = 0
C
C             Book JDRAW structure for view banks
C
         CALL MZBOOK(IXCONS,JDRAW,JDRAW,1,'DRAW',0,0,0,3,0)
      ENDIF
 
C
      IF(JHEAD.GT.0)THEN
         IDRUN=IQ(JHEAD+1)
         IDEVT=IQ(JHEAD+2)
      ENDIF
      IF(JRUNG.GT.0) THEN
*
*             Here we deal with version numbers If JRUNG has been read in,
*             then save the version numbers of the new JRUNG and restore
*             the current version number for KINE, HITS and DIGI
         DO 300 J=1, NLINK
            IF(IVERSI(J).GT.0) THEN
               NKEY = ABS(LINK(J))
               IF(KNAMES(NKEY).EQ.'RUNG') THEN
                  DO 280 I=1,8
                     IDNEW(I) = IQ(JRUNG+10+I)
                     VERNEW(I) = Q(JRUNG+20+I)
  280             CONTINUE
*
*             And we put back the old version numbers because,
*             in principle, KINE, HITS and DIGI have not be read in
                  DO 290 I=3,8
                     IQ(JRUNG+10+I) = IDOLD(I)
                     Q(JRUNG+20+I) = VEROLD(I)
  290             CONTINUE
               ENDIF
            ENDIF
  300    CONTINUE
*
*            And here we do it again for KINE, HITS and DIGI
         DO 310 J=1, NLINK
            IF(IVERSI(J).GT.0) THEN
               NKEY = ABS(LINK(J))
               IF(KNAMES(NKEY).EQ.'KINE') THEN
                  IF(IDNEW(3).GT.0) THEN
                     IQ(JRUNG+13) = IDNEW(3)
                     IQ(JRUNG+14) = IDNEW(4)
                     Q(JRUNG+23) = VERNEW(3)
                     Q(JRUNG+24) = VERNEW(4)
                  ENDIF
               ELSEIF(KNAMES(NKEY).EQ.'HITS') THEN
                  IF(IDNEW(5).GT.0) THEN
                     IQ(JRUNG+15) = IDNEW(5)
                     IQ(JRUNG+16) = IDNEW(6)
                     Q(JRUNG+25) = VERNEW(5)
                     Q(JRUNG+26) = VERNEW(6)
                  ENDIF
               ELSEIF(KNAMES(NKEY).EQ.'DIGI') THEN
                  IF(IDNEW(7).GT.0) THEN
                     IQ(JRUNG+17) = IDNEW(7)
                     IQ(JRUNG+18) = IDNEW(8)
                     Q(JRUNG+27) = VERNEW(7)
                     Q(JRUNG+28) = VERNEW(8)
                  ENDIF
               ELSEIF(KNAMES(NKEY).EQ.'MATE'.OR. KNAMES(NKEY) .EQ.'TMED'
     +         ) THEN
                  IF(VERNEW(1).NE.0) THEN
*                We know which version number we are reading
                     IF(VERNEW(1).LT.GVERSN) THEN
                        WRITE(CHMAIL,10900) KNAMES(NKEY),VERNEW(1),
     +                  GVERSN
                        IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
                        WRITE(CHMAIL,11000)
                        IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  310    CONTINUE
      ENDIF
  320 CONTINUE
*
10000 FORMAT(' *** GFIN *** Key ',A4,' ignored for initialization')
10100 FORMAT(' *** GFIN *** Key ',A4,' ignored for kinematics')
10200 FORMAT(' *** GFIN *** Key ',A4,' ignored for trigger')
10300 FORMAT(' *** GFIN *** Unknown key ',A4)
10400 FORMAT(' *** GFIN *** No valid key given')
10500 FORMAT(' *** GFIN *** Data structure ',A4,' version ',I10,
     +       ' successfully read in ')
10600 FORMAT(' *** GFIN *** Data structure ',A4,' was not found')
10700 FORMAT(' *** GFIN *** Data structure ',A4,' version ',I10,
     +       ' was not found')
10800 FORMAT(' *** GFIN *** Nothing found to read !')
10900 FORMAT(' *** GFIN *** ',A4,' data structure ',
     +       'version ',F6.4,' current version is ',F6.4)
11000 FORMAT('              Please call subroutine GPHYSI before ',
     +       'tracking')
11100 FORMAT(' *** GFIN ***  Illegal number of links ',I10)
  999 END

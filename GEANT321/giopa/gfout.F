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
      SUBROUTINE GFOUT(LUN,CHOBJ,NKEYS,IDVERS,CHOPT,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to write GEANT object(s) into the FZ file        *
C.    *       The data structures in memory are written on disk        *
C.    *                                                                *
C.    *       LUN    Logical unit                                      *
C.    *                                                                *
C.    *       CHOBJ  The type of data structure to be written:         *
C.    *              MATE material                                     *
C.    *              TMED tracking medium                              *
C.    *              VOLU volumes                                      *
C.    *              ROTM rotation matrix                              *
C.    *              SETS detector set                                 *
C.    *              PART particle                                     *
C.    *              SCAN geometry                                     *
C.    *              INIT all above                                    *
C.    *              KINE this keyword will trigger the write of       *
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
C.    *              TRIG this keyword will trigger the write of       *
C.    *                   DIGI, HEAD, HITS, KINE, VERT abd JXYZ        *
C.    *                   unless the 'S' flag is set                   *
C.    *                                                                *
C.    *       NKEYS  number of keys in vector CHOBJ                    *
C.    *                                                                *
C.    *       IDVERS version of the data structure to be written out   *
C.    *                                                                *
C.    *       CHOPT  List of options                                   *
C.    *                   'I'      write only initialisation data      *
C.    *                            structures                          *
C.    *                   'K'      write only KINE and TRIG data       *
C.    *                            structures                          *
C.    *                   'T'      write only DIGI, HEAD, HITS, KINE,  *
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
C.    *    The FZ data base can be read in via GOPEN/GFIN              *
C.    *                                                                *
C.    *                                                                *
C.    *      Example.                                                  *
C.    *                                                                *
C.    *      CALL GOPEN(1,'O',1024,IER)                                *
C.    *      CALL GFOUT (1,'VOLU',1,0,' ',IER)                         *
C.    *      CALL GFOUT (1,'MATE',1,0,' ',IER)                         *
C.    *      CALL GFOUT (1,'TMED',1,0,' ',IER)                         *
C.    *      CALL GFOUT (1,'ROTM',1,0,' ',IER)                         *
C.    *      CALL GFOUT (1,'PART',1,0,' ',IER)                         *
C.    *      CALL GFOUT (1,'SCAN',1,0,' ',IER)                         *
C.    *      CALL GFOUT (1,'SETS',1,0,' ',IER)                         *
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
      DIMENSION IUHEAD(2),ITRAN(23),JTRAN(23)
      EQUIVALENCE (JNAMES(1),JDIGI)
      CHARACTER*4 KNAMES(NMKEY),CHOBJ(*)
      CHARACTER*(*) CHOPT
      DATA KNAMES/'DIGI','DRAW','HEAD','HITS','KINE','MATE','PART',
     +     'ROTM','RUNG','SETS','STAK','STAT','TMED','NULL','VERT',
     +     'VOLU','JXYZ','NULL','NULL','NULL','SCAN','NULL'/
      DATA ITRAN/7,6,13,16,8,10,2,9,8*0,3,15,5,17,4,1,21/
      DATA JTRAN/22,7,17,21,19,2,1,5,8,6,2*0,3,0,18,4,20,3*0,23,2*0/
      DATA IXD/2,1,2,2,2,8*1,2,2,1,2,3*0,1,0/
      DATA LINIT/2,6,7,8,9,10,13,16,21/
      DATA LKINE/5,15/
      DATA LTRIG/1,3,4,5,15,17/
C.
C.    ------------------------------------------------------------------
C.
      IQUEST(1)=0
      IER=0
      LDIV(1)  =IXCONS
      LDIV(2)  =IXDIV
*
      IOPTI=INDEX(CHOPT,'i')+INDEX(CHOPT,'I')
      IOPTT=INDEX(CHOPT,'t')+INDEX(CHOPT,'T')
      IOPTK=INDEX(CHOPT,'k')+INDEX(CHOPT,'K')
      IOPTS=INDEX(CHOPT,'s')+INDEX(CHOPT,'S')
      IOPTQ=INDEX(CHOPT,'q')+INDEX(CHOPT,'Q')
*
      NLINK=0
      DO 90  JKEY=1,NKEYS
         IF(IOPTS.EQ.0) THEN
         IF(CHOBJ(JKEY).EQ.'INIT') THEN
            DO 20 J=1, NLINIT
               DO 10  MLINK=1,NLINK
                  IF(LINK(MLINK).EQ.LINIT(J)) GO TO 20
   10          CONTINUE
               NLINK=NLINK+1
               LINK(NLINK)=LINIT(J)
   20       CONTINUE
            GO TO 90
         ELSEIF(CHOBJ(JKEY).EQ.'TRIG') THEN
            DO 40 J=1, NLTRIG
               DO 30  MLINK=1,NLINK
                  IF(LINK(MLINK).EQ.LTRIG(J)) GO TO 40
   30          CONTINUE
               NLINK=NLINK+1
               LINK(NLINK)=LTRIG(J)
   40       CONTINUE
            GO TO 90
         ELSEIF(CHOBJ(JKEY).EQ.'KINE') THEN
            DO 60 J=1, NLKINE
               DO 50  MLINK=1,NLINK
                  IF(LINK(MLINK).EQ.LKINE(J)) GO TO 60
   50          CONTINUE
               NLINK=NLINK+1
               LINK(NLINK)=LKINE(J)
   60       CONTINUE
            GO TO 90
         ENDIF
         ENDIF
            DO 80 J=1,NMKEY
               IF(CHOBJ(JKEY).EQ.KNAMES(J)) THEN
                  DO 70 MLINK=1,NLINK
                     IF(LINK(MLINK).EQ.J) GO TO 90
   70             CONTINUE
                  NLINK=NLINK+1
                  LINK(NLINK)=J
                  GO TO 90
               ENDIF
   80       CONTINUE
            WRITE(CHMAIL,10300) CHOBJ(JKEY)
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
   90 CONTINUE
*
      IF(IOPTI.GT.0) THEN
         DO 110 J=1, NLINK
            DO 100 K=1, NLINIT
               IF(LINK(J).EQ.LINIT(K)) GO TO 110
  100       CONTINUE
            WRITE(CHMAIL,10000) KNAMES(LINK(J))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            LINK(J)=0
  110    CONTINUE
      ELSEIF(IOPTK.GT.0) THEN
         DO 130 J=1, NLINK
            DO 120 K=1, NLKINE
               IF(LINK(J).EQ.LKINE(K)) GO TO 130
  120       CONTINUE
            WRITE(CHMAIL,10100) KNAMES(LINK(J))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            LINK(J)=0
  130    CONTINUE
      ELSEIF(IOPTT.GT.0) THEN
         DO 150 J=1, NLINK
            DO 140 K=1, NLTRIG
               IF(LINK(J).EQ.LTRIG(K)) GO TO 150
  140       CONTINUE
            WRITE(CHMAIL,10200) KNAMES(LINK(J))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            LINK(J)=0
  150    CONTINUE
      ENDIF
*
      IOFF=0
      DO 160 J=1, NLINK
         IF(LINK(J).EQ.0) THEN
           IOFF=IOFF-1
         ELSE
           LINK(J+IOFF)=LINK(J)
         ENDIF
  160 CONTINUE
      NLINK=NLINK+IOFF
      IF(IOPTI+IOPTK+IOPTT.EQ.0) THEN
*
*        We have to choose which event header to write, JRUNG or JHEAD
*        If the banks list contains banks which depends on both headers,
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
  163      CONTINUE
  168    CONTINUE
  169    CONTINUE
      ENDIF
*
      IF(NLINK.LE.0) THEN
         WRITE(CHMAIL,10400)
         IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         IER=-1
         GOTO 999
      ENDIF
*
      NWOUT=0
      IOFW =0
      DO 170 J=1,NLINK
         IVERSI(J)=0
         IRESUL(J)=0
         NKEY=LINK(J)
         IF(NKEY.EQ.3.OR.NKEY.EQ.9) THEN
            IOFW=1
            NPOS=J
         ENDIF
         LINK(J)=-ABS(LINK(J))
         IF(NKEY.LE.20)THEN
            IF(JNAMES(NKEY).NE.0) THEN
               LINK(J)=ABS(LINK(J))
               NWOUT=NWOUT+1
            ENDIF
         ELSE
            NKL=NKEY-20
            IF(ISLINK(NKL).NE.0) THEN
               LINK(J)=ABS(LINK(J))
               NWOUT=NWOUT+1
            ENDIF
         ENDIF
  170 CONTINUE
*
*               Write next start of event data structure
      IUHEAD(1)=IDVERS
      IUHEAD(2)=NWOUT-IOFW
      NUH=2
      IF(IOPTI.NE.0) THEN
         CALL FZOUT(LUN,IXCONS,JRUNG,1,'L',2,NUH,IUHEAD)
      ELSEIF(IOPTT+IOPTK.NE.0) THEN
         CALL FZOUT(LUN,IXDIV,JHEAD,1,'L',2,NUH,IUHEAD)
      ENDIF
      IF(IQUEST(1).EQ.0) THEN
         IVERSI(NPOS)=IDVERS
         IRESUL(NPOS)=1
      ELSE
         WRITE(CHMAIL,10500) KNAMES(LINK(NPOS))
      ENDIF
*
      DO 180 IK=1,NLINK
*
*              Write selected data structures
         NKEY=LINK(IK)
         IF(NKEY.GT.0) THEN
            IF(NKEY.EQ.9) THEN
               GOTO 180
            ELSEIF(NKEY.EQ.3) THEN
               GOTO 180
            ELSEIF(NKEY.EQ.1) THEN
               CALL GRLEAS(JDIGI)
            ELSEIF(NKEY.EQ.4) THEN
               CALL GRLEAS(JHITS)
            ENDIF
            IDIV=LDIV(IXD(NKEY))
            JKEY=JTRAN(NKEY)
            IF(NKEY.LE.20)THEN
               CALL FZOUT(LUN,IDIV,JNAMES(NKEY),0,'L',2,1,JKEY)
            ELSE
               NKL=NKEY-20
               CALL FZOUT(LUN,IDIV,ISLINK(NKL),0,'L',2,1,JKEY)
            ENDIF
            IF(IQUEST(1).EQ.0) THEN
               IVERSI(IK)=IDVERS
               IRESUL(IK)=1
            ELSE
               WRITE(CHMAIL,10500) KNAMES(NKEY)
            ENDIF
         ENDIF
  180 CONTINUE
*
      NOUT=0
      DO 190 I=1,NLINK
         IF(IRESUL(I).EQ.1) THEN
            WRITE(CHMAIL,10600) KNAMES(ABS(LINK(I))),IVERSI(I)
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
            NOUT=NOUT+1
         ELSE
            WRITE(CHMAIL,10700) KNAMES(ABS(LINK(I)))
            IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         ENDIF
  190 CONTINUE
*
      IF(NOUT.LE.0) THEN
         WRITE(CHMAIL,10800)
         IF(IOPTQ.EQ.0) CALL GMAIL(0,0)
         IER=-1
      ELSEIF(NOUT.LT.NLINK) THEN
         IER=NOUT
      ENDIF
*
10000 FORMAT(' *** GFOUT *** Key ',A4,' ignored for initialization')
10100 FORMAT(' *** GFOUT *** Key ',A4,' ignored for kinematics')
10200 FORMAT(' *** GFOUT *** Key ',A4,' ignored for trigger')
10300 FORMAT(' *** GFOUT *** Unknown key ',A4)
10400 FORMAT(' *** GFOUT *** No valid key given')
10500 FORMAT(' *** GFOUT *** Problems writing data structure ',A4)
10600 FORMAT(' *** GFOUT *** Data structure ',A4,' version ',I10,
     +       ' successfully written out')
10700 FORMAT(' *** GFOUT *** Data structure ',A4,' not found')
10800 FORMAT(' *** GFOUT *** Nothing written out !')
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
*
      SUBROUTINE GDPLST(JVSIM,NAME,LEVMAX,KXXX)
*
***************************************************************
*                                                             *
*  Routine name: GDPLST                     Date: 19-2-92     *
*                                                             *
*  Function: It performs the drawing of the tree under the    *
*            user control. Clicking the left button of the    *
*            mouse when the cursor is inside the name's box   *
*            will display the specification of such a volume, *
*            clicking on the red arrows will display the tree *
*            below it, clicking on the green arrows will dis- *
*            play the tree starting from an upper level.      *
*            The number of levels displayed in both cases     *
*            depends on the arrow which has been clicked.     *
*                                                             *
*  Author: S. Giani   ******                                  *
*                                                             *
***************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/pawc.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gcmutr.inc"
#include "geant321/gcspee.inc"
#include "geant321/gchil2.inc"
#include "geant321/gccurs.inc"
#include "geant321/gcursb.inc"
*
      CHARACTER*4 NAME,SAVNAM,NAME1,MOTH,SAVSAN,NNNN
      CHARACTER*1 GCHAR
*
      IF(NNPAR.EQ.3.AND.KXXX.EQ.0)GOTO 110
      JSIM=0
      IDINAM=2
      IF(JAASS1.NE.0) CALL MZDROP(IXSTOR, JAASS1, ' ')
      IF(JAASS2.NE.0) CALL MZDROP(IXSTOR, JAASS2, ' ')
      IF(JAASS3.NE.0) CALL MZDROP(IXSTOR, JAASS3, ' ')
      IF(JAASS4.NE.0) CALL MZDROP(IXSTOR, JAASS4, ' ')
      IF(JTICKS.NE.0) CALL MZDROP(IXSTOR, JTICKS, ' ')
      IF(JMYLLS.NE.0) CALL MZDROP(IXSTOR, JMYLLS, ' ')
      IF(JMYMOT.NE.0) CALL MZDROP(IXSTOR, JMYMOT, ' ')
      NEE=NUMND1+10
      NEES=NEE*7
      CALL MZNEED(IXDIV,NEES,'G')
      CALL MZBOOK(IXDIV,JAASS1,JAASS1,1,'TREX',0,0,NEE,3,-1)
      CALL MZBOOK(IXDIV,JAASS2,JAASS2,1,'TREY',0,0,NEE,3,-1)
      CALL MZBOOK(IXDIV,JAASS3,JAASS3,1,'NAMS',0,0,NEE,2,-1)
      CALL MZBOOK(IXDIV,JAASS4,JAASS4,1,'SCAS',0,0,NEE,2,-1)
      CALL MZBOOK(IXDIV,JTICKS,JTICKS,1,'TICK',0,0,NEE,2,-1)
      CALL MZBOOK(IXDIV,JMYLLS,JMYLLS,1,'MYLS',0,0,NEE,2,-1)
      CALL MZBOOK(IXDIV,JMYMOT,JMYMOT,1,'MYMO',0,0,NEE,2,-1)
      CALL IGRNG(YPLT,XPLT)
      SAVNAM=NAME
      SAVSAN=SAVNAM
      LSAVLE=LEVMAX
      LSAVSA=LSAVLE
      YPLTNE=YPLT
      XPLTNE=XPLT
      NUMNDS=NUMND1
      DO 10  J=1,NUMND1
         Q(JAASS1+J)=Q(JULEV+J)
         Q(JAASS2+J)=Q(JVLEV+J)
         IQ(JAASS3+J)=IQ(JNAM1+J)
         IQ(JAASS4+J)=IQ(JSCA1+J)
         IQ(JTICKS+J)=IQ(JTICK+J)
         IQ(JMYLLS+J)=IQ(JMYLL+J)
         IQ(JMYMOT+J)=IQ(JMOT1+J)
   10 CONTINUE
      IPLAC=LEVMAX
   20 CONTINUE
      IPICK=0
      IF(IDINAM.EQ.2)THEN
        CALL IGRNG(YPLT,XPLT)
      ELSE
        CALL IGRNG(YPLTNE,XPLTNE)
      ENDIF
      IF(KXXX.EQ.1)THEN
        KACKT=1
        KCHAR=1
        IPICK=1
         CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
         IF (IVO.LE.0) GO TO 110
        GOTO 50
      ENDIF
      CALL IRQLC(1,1,KCHAR,NT,V0,U0)
      IF (KCHAR.EQ.0) GO TO 50
      GCHAR=CHAR(KCHAR)
      IF(IDINAM.EQ.2)THEN
        NUMNUM=NUMND1
      ELSE
        NUMNUM=NUMNDS
      ENDIF
      DO 40  J=1,NUMNUM
       IF(IDINAM.EQ.2)THEN
         IF (IQ(JSCA1+J).NE.0) THEN
            IF(U0.LT.(Q(JULEV+J)-SIZD2).OR.U0.GT.(Q(JULEV+J)+SIZD2))
     +      GO TO 30
            IF(V0.LT.(Q(JVLEV+J)-(SIZD2*2.*FACHV)-HALF).OR.V0.GT.(
     +      Q(JVLEV+J)+(SIZD2*2.*FACHV)+HALF)) GO TO 30
            IPICK=J
            KACKT=0
            IF(V0.LT.(Q(JVLEV+J)-(SIZD2*2.*FACHV)))THEN
               KACKT=-1
               PLAC=((Q(JVLEV+J)-(SIZD2*2.*FACHV))-V0)*(IQ(JMYLL+J)
     +         -1)/HALF
               IPLACS=PLAC+1
            ENDIF
            IF(V0.GT.(Q(JVLEV+J)+(SIZD2*2.*FACHV)))THEN
               KACKT=1
               PLAC=(V0-(Q(JVLEV+J)+(SIZD2*2.*FACHV)))*IQ(JTICK+J)/
     +         HALF
               IPLAC=PLAC
*                 IF((PLAC-IPLAC).GT..5)IPLAC=IPLAC+1
               IPLAC=IPLAC+2
            ENDIF
            IF (KCHAR.EQ.1.OR.KCHAR.EQ.2) GO TO 50
            IF(GCHAR.EQ.' ')GO TO 50
            CALL UHTOC(IQ(JNAM1+J),4,NAME1,4)
            IF (GCHAR.EQ.'2') THEN
               CALL GSATT(NAME1,'SEEN',-2)
            ELSE IF (GCHAR.EQ.'-') THEN
               CALL GSATT(NAME1,'SEEN',-1)
            ELSE IF (GCHAR.EQ.'0') THEN
               CALL GSATT(NAME1,'SEEN',0)
            ELSE IF (GCHAR.EQ.'1') THEN
               CALL GSATT(NAME1,'SEEN',1)
            ENDIF
            GO TO 20
         ENDIF
       ENDIF
   30    CONTINUE
   40 CONTINUE
   50 CONTINUE
      IF(KCHAR.EQ.0)THEN
       IF(IDINAM.EQ.2)THEN
         IF(JVSIM.EQ.2.AND.JSIM.EQ.1)CALL ICLWK(2)
         GOTO 110
       ENDIF
      ELSE
         IF (IPICK.NE.0) THEN
          IF(JSIM.EQ.0.AND.JVSIM.EQ.2.AND.KACKT.EQ.0)THEN
            JSIM=1
            CALL IOPWK(2,1,9)
          ENDIF
           IF(IDINAM.EQ.1)THEN
            CALL UHTOC(IQ(JAASS3+IPICK),4,NAME,4)
            CALL IDAWK(1)
            CALL IACWK(JVSIM)
           ELSEIF(IDINAM.EQ.2.AND.KXXX.NE.1)THEN
            CALL UHTOC(IQ(JNAM1+IPICK),4,NAME,4)
           ENDIF
            CALL ISELNT(1)
*            CALL ICLRWK(JVSIM,1)
            IF(KACKT.NE.0)CALL ICLRWK(0,0)
            IF(KACKT.EQ.0)THEN
               CALL IDAWK(1)
               CALL IACWK(JVSIM)
               CALL IGRNG(SAVPLX,SAVPLY)
               CALL ICLRWK(0,0)
               CALL GDSPEC(NAME)
               CALL IDAWK(JVSIM)
               CALL IACWK(1)
               IF(IDINAM.EQ.2.AND.JVSIM.EQ.1)THEN
                     CALL IRQLC(1,1,KCHAR,NT,V0,U0)
                     CALL ICLRWK(0,0)
                     CALL ISELNT(1)
                     INTFLA=1
                     CALL GDTREE(SAVSAN,0,110)
                     INTFLA=-1
                     CALL GDTREE(SAVSAN,LSAVSA,111)
               ENDIF
            ELSEIF(KACKT.EQ.1)THEN
               IF(KXXX.EQ.1)KXXX=0
               NUMNDA=NUMND1
               NUMND1=NUMND2
               NNNN=NAME
               IADDI=0
               DO 81  JI=1,NUMND2
                  IQ(JNAM1+JI)=IQ(JFINAM+JI)
                  IQ(JSCA1+JI)=IQ(JFISCA+JI)
                  IQ(JMOT1+JI)=IQ(JFIMOT+JI)
   81          CONTINUE
               IF(NAME.EQ.MOMO)GOTO 41
               DO 91  II=1,16
                  CALL GDTR8(NAME,MOTH,IONL)
                  IF(MOTH.EQ.MOMO)THEN
                   IADDI=II
                   GOTO 41
                  ENDIF
                  NAME=MOTH
   91          CONTINUE
   41          CONTINUE
               NUMND1=NUMNDA
               NAME=NNNN
               INTFLA=1
               CALL GDTREE(NAME,0,110)
               INTFLA=-1
               IF(NNPAR.NE.3)THEN
                CALL GDTREE(NAME,IPLAC,111)
               ELSE
                CALL GDTREE(NAME,IPLAC,IISELT)
                GOTO 110
               ENDIF
               IF(IDINAM.EQ.2.OR.KCHAR.EQ.2)THEN
                     SAVSAN=NAME
                     LSAVSA=IPLAC
               ENDIF
            ELSEIF(KACKT.EQ.-1)THEN
               NUMNDA=NUMND1
               NUMND1=NUMND2
               IADDI=0
               DO 80  JII=1,NUMND2
                  IQ(JNAM1+JII)=IQ(JFINAM+JII)
                  IQ(JSCA1+JII)=IQ(JFISCA+JII)
                  IQ(JMOT1+JII)=IQ(JFIMOT+JII)
   80          CONTINUE
               DO 90  II=1,IPLACS
                  CALL GDTR8(NAME,MOTH,IONL)
                  NAME=MOTH
   90          CONTINUE
               NNNN=NAME
               IF(NAME.EQ.MOMO)GOTO 42
               DO 92  II=1,16
                  CALL GDTR8(NAME,MOTH,IONL)
                  IF(MOTH.EQ.MOMO)THEN
                   IADDI=II
                   GOTO 42
                  ENDIF
                  NAME=MOTH
   92          CONTINUE
   42          CONTINUE
               NAME=NNNN
               NUMND1=NUMNDA
               INTFLA=1
               CALL GDTREE(NAME,0,110)
               INTFLA=-1
               CALL GDTREE(NAME,IPLAC,111)
               IF(IDINAM.EQ.2.OR.KCHAR.EQ.2)THEN
                     SAVSAN=NAME
                     LSAVSA=IPLAC
               ENDIF
            ENDIF
            CALL ISELNT(1)
            GO TO 20
         ELSE
            GOTO 20
         ENDIF
      ENDIF
*
 110  CONTINUE
***************      CALL IGRNG(SAVPLX,SAVPLY)
      IF(JAASS1.NE.0) CALL MZDROP(IXSTOR, JAASS1, ' ')
      IF(JAASS2.NE.0) CALL MZDROP(IXSTOR, JAASS2, ' ')
      IF(JAASS3.NE.0) CALL MZDROP(IXSTOR, JAASS3, ' ')
      IF(JAASS4.NE.0) CALL MZDROP(IXSTOR, JAASS4, ' ')
      IF(JTICKS.NE.0) CALL MZDROP(IXSTOR, JTICKS, ' ')
      IF(JMYLLS.NE.0) CALL MZDROP(IXSTOR, JMYLLS, ' ')
      IF(JMYMOT.NE.0) CALL MZDROP(IXSTOR, JMYMOT, ' ')
      IF(JMYLL.NE.0) CALL MZDROP(IXSTOR, JMYLL, ' ')
      IF(JTICK.NE.0) CALL MZDROP(IXSTOR, JTICK, ' ')
               IF(JFIMOT.NE.0) CALL MZDROP(IXSTOR, JFIMOT, ' ')
               IF(JFISCA.NE.0) CALL MZDROP(IXSTOR, JFISCA, ' ')
               IF(JFINAM.NE.0) CALL MZDROP(IXSTOR, JFINAM, ' ')
      LARETT(1)=0
  999 END

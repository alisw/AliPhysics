*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTREE(KNAME,LEVMAX,ISEL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw the tree of geometric data structure starting       *
C.    *       from object KNAME, for LEVMAX depth levels (0=all)       *
C.    *                                                                *
C.    *       If ISEL=      0   then draw only node name;              *
C.    *       if ISEL=xxxxxx1   then add multiplicity;                 *
C.    *       if ISEL=xxxxx1x   then add 'ONLY' information;           *
C.    *       if ISEL=xxxx1xx   then add 'DET ' information;           *
C.    *       if ISEL=xxx1xxx   then add 'SEEN' information;           *
C.    *       if ISEL=xx1xxxx   then add little picture of volume;     *
C.    *       if ISEL=x1xxxxx   then add graphics cursor and           *
C.    *                         returns in IPICK the node picked;      *
C.    *                                                                *
C.    *       For very big trees, the attribute SEEN -3 can be applied *
C.    *       to any node in the level that has too many different     *
C.    *       nodes (each with a different volume name) : the routine  *
C.    *       GDTREE will draw only the first, the last, and one       *
C.    *       dummy node in the middle that tells how many nodes       *
C.    *       does it stand for.                                       *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Authors : P.Zanarini ; S.Giani  *********                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gchil2.inc"
#include "geant321/gccurs.inc"
#include "geant321/gcursb.inc"
*
      CHARACTER*4 KNAME,NAME,NAME1,MOTH
      DIMENSION U(4),V(4),UUU(15),VVV(15),SVVX(3),SVVY(3)
      SAVE SAVTHE,SAVPHI,SAVU0,SAVV0,SAVSCU,SAVSCV,SVGZUA,SVGZVA
      SAVE SVGZUB,SVGZVB,SVGZUC,SVGZVC
      COMMON/NPILEV/NPILE
 
C.
C.    ------------------------------------------------------------------
C.
C             Is NAME an existing volume ?
C
*      CALL IGSET('SYNC',1.)
      IF (MOMO.EQ.' ') THEN
         CALL UHTOC(IQ(JVOLUM+1),4,NAME,4)
         MOMO=NAME
      ENDIF
      IF (KNAME.EQ.' ') THEN
         CALL UHTOC(IQ(JVOLUM+1),4,NAME,4)
         KNAME=NAME
      ELSE
         CALL GLOOK(KNAME,IQ(JVOLUM+1),NVOLUM,IVO)
         IF (IVO.LE.0) GO TO 180
         NAME=KNAME
      ENDIF
C
C             Build tree structure using view bank 11
C
      IVTREE=NKVIEW+1
      CALL UCTOH (NAME, IROOT, 4, 4)
      CALL GDTR0 (IVTREE, IROOT, LEVMAX, IER)
      IF (IER.NE.0) GO TO 170
C
      IPICK=0
C
      ISEL5=ISEL/100000
      ISEL5=ISEL5-ISEL5/2*2
      ISEL4=ISEL/10000
      ISEL4=ISEL4-ISEL4/2*2
      ISEL3=ISEL/1000
      ISEL3=ISEL3-ISEL3/2*2
      ISEL2=ISEL/100
      ISEL2=ISEL2-ISEL2/2*2
      ISEL1=ISEL/10
      ISEL1=ISEL1-ISEL1/2*2
      ISEL0=ISEL/1
      ISEL0=ISEL0-ISEL0/2*2
C
C             Save GDRAW calling parameters
C             and ZOOM internal parameters
C
      IDRAW=0
      IF (ISEL4.NE.0) THEN
         IDRAW=1
         SAVTHE=GTHETA
         SAVPHI=GPHI
         SAVU0=GU0
         SAVV0=GV0
         SAVSCU=GSCU
         SAVSCV=GSCV
         SVGZUA=GZUA
         SVGZVA=GZVA
         SVGZUB=GZUB
         SVGZVB=GZVB
         SVGZUC=GZUC
         SVGZVC=GZVC
         GZUA=1
         GZVA=1
         GZUB=0
         GZVB=0
         GZUC=0
         GZVC=0
      ENDIF
C
C             Fill arrays Q(JULEV) and Q(JVLEV)
C
      DO 10 I=1,NUMND1
         IQ(JSCA1+I)=0
   10 CONTINUE
C
      MAXV=LEVMAX
      LEVVER=1
      LEVHOR=0
      MLEVV=1
      CALL GDTR10(1)
      MLEVH=LEVHOR
C
      PLTVER=25.
      SIZE=PLTVER/(4.*MLEVV)
***      SIZE=PLTVER/(4.*MLEVV)
C
C             Compute user coordinates boundaries of tree picture
C
      SIZ2=SIZE*2.
***      SIZ4=SIZE*4.
      IF(MLEVH.GT.MLEVV)THEN
         SIZ4=(MLEVH*SIZ2)/MLEVV
         FACHV=(MLEVH/2.)/MLEVV
         IF(FACHV.LT.2.)FACHV=1.
      ELSE
         SIZ4=SIZE*4.
         FACHV=1.
      ENDIF
      SIZD2=SIZE/2.
      XPLT=MLEVH*SIZ2
      YPLT=MLEVV*SIZ4
C
C             Save current ranges and store new ones
C
      SAVPLX=PLTRNX
      SAVPLY=PLTRNY
C
      PLTRNX=YPLT
      PLTRNY=XPLT
C
      CALL IGRNG(YPLT,XPLT)
*      CALL IGRNG(XPLT,YPLT)
C
C             Draw nodes
C
      IDUP=1
      IONL=1
      IDET=0
      ISEEN=1
*      INTSPI=0
      DO 60 J=1,NUMND1
         IF (IQ(JSCA1+J).NE.0) THEN
            CALL UHTOC(IQ(JNAM1+J),4,NAME,4)
            IF(INTFLA.EQ.10)THEN
               IADDI=0
*               INTSPI=1
               INTFLA=1
               CALL MZLINT(IXDIV,'/GCHIL2/',LARETT,JMYMOT,LARETT)
               LARETT(1)=1
               IF(JTICK.NE.0) CALL MZDROP(IXSTOR, JTICK, ' ')
               IF(JMYLL.NE.0) CALL MZDROP(IXSTOR, JMYLL, ' ')
               NEE=NUMND1+10
               NEES=NEE*2
               CALL MZNEED(IXDIV,NEES,'G')
               CALL MZBOOK(IXDIV,JTICK,JTICK,1,'TTT',0,0,NEE,2,-1)
               CALL MZBOOK(IXDIV,JMYLL,JMYLL,1,'MMM',0,0,NEE,2,-1)
               IF(JFIMOT.NE.0) CALL MZDROP(IXSTOR, JFIMOT, ' ')
               IF(JFISCA.NE.0) CALL MZDROP(IXSTOR, JFISCA, ' ')
               IF(JFINAM.NE.0) CALL MZDROP(IXSTOR, JFINAM, ' ')
               NEES=NEE*3
               CALL MZNEED(IXDIV,NEES,'G')
               CALL MZBOOK(IXDIV,JFIMOT,JFIMOT,1,'FIMO',0,0,NEE,2,-1)
               CALL MZBOOK(IXDIV,JFISCA,JFISCA,1,'FISC',0,0,NEE,2,-1)
               CALL MZBOOK(IXDIV,JFINAM,JFINAM,1,'FINA',0,0,NEE,2,-1)
               NUMND2=NUMND1
            ENDIF
            IF(INTFLA.EQ.1)THEN
               MLETMP=MLEVV
               LEVTMP=LEVVER
               MLEVV=1
               LEVVER=1
               CALL GDTR10(J)
               IQ(JTICK+J)=MLEVV-1
               MLEVV=MLETMP
               LEVVER=LEVTMP
               IF(J.NE.1)THEN
*               IF(INTSPI.NE.1)THEN
*               NUMNDA=NUMND1
*               NUMND1=NUMNDS
*               DO 89  JI=1,NUMNDS
*                  IQ(JNAM1+JI)=IQ(JAASS3+JI)
*                  IQ(JSCA1+JI)=IQ(JAASS4+JI)
*                  IQ(JMOT1+JI)=IQ(JMYMOT+JI)
*   89          CONTINUE
*               ENDIF
                  DO 20  LL=1,16
                     CALL GDTR8(NAME,MOTH,IONL)
                     IF(MOTH.EQ.KNAME)THEN
                        IQ(JMYLL+J)=LL+1+IADDI
                        GOTO 30
                     ENDIF
                     NAME=MOTH
   20             CONTINUE
   30             CONTINUE
*                  IF(INTSPI.NE.1)NUMND1=NUMNDA
               ELSE
                  IQ(JMYLL+J)=1+IADDI
               ENDIF
            ELSE
               Q(JULEV+J) = (Q(JULEV+J)-1.)*SIZ2 + SIZE
***            Q(JVLEV+J) = YPLT - (Q(JVLEV+J)-1.)*SIZ4 - SIZ2
               Q(JVLEV+J) = (Q(JVLEV+J)-1.)*SIZ4 + SIZ2
               IF (ISEL0.NE.0) IDUP=IQ(JDUP1+J)
               IF (ISEL1.NE.0) CALL GDTR8(NAME,MOTH,IONL)
               IF (ISEL2.NE.0) THEN
                  IDET=0
                  KVAL=0
                  CALL GFATT(NAME,'DET ',KVAL)
                  IF (KVAL.GT.0) IDET=1
               ENDIF
               IF (ISEL3.NE.0) THEN
                  KVAL=1
                  CALL GFATT(NAME,'SEEN',KVAL)
                  IF (KVAL.EQ.1.OR.KVAL.EQ.-2) THEN
                     ISEEN=1
                  ELSE
                     ISEEN=0
                  ENDIF
                  IF (J.EQ.1) GO TO 50
                  JM=J
   40             CONTINUE
                  JM=IQ(JMOT1+JM)
                  KVALM=1
                  CALL UHTOC(IQ(JNAM1+JM),4,NAME1,4)
                  CALL GFATT(NAME1,'SEEN',KVALM)
                  IF (KVALM.LT.0) THEN
                     ISEEN=0
                     GO TO 50
                  ENDIF
                  IF (JM.NE.1) GO TO 40
   50             CONTINUE
               ENDIF
*****            CALL GDPRTR(NAME,Q(JVLEV+J),Q(JULEV+J),SIZE,IDUP,IDRAW,
*****     +      IONL, IDET,ISEEN)
               NPILE=IQ(JMYLL+J)
               CALL GDPRTR(NAME,Q(JULEV+J),Q(JVLEV+J),SIZE,FACHV,IDUP,
     +         IDRAW, IONL, IDET,ISEEN)
            ENDIF
         ENDIF
   60 CONTINUE
C
C             Draw links
C
      IF(INTFLA.EQ.1)GOTO 160
      LINCOL=1
      CALL MVBITS(LINCOL,0,8,LINATT,16)
      IF (IQ(JMOT1+1).NE.0) THEN
         U(2)=Q(JULEV+1)
         U(1)=U(2)
         V(2)=Q(JVLEV+1)-(SIZD2*2.*FACHV)
         V(1)=V(2)-(SIZE*2.*FACHV)
         CALL GDRAWV(V,U,2)
***         CALL GDRAWV(U,V,2)
      ENDIF
C
      IGREEN=1
      DO 150 J=1,NUMND1
         IF (IQ(JSCA1+J).EQ.1) THEN
            JX=IQ(JXON1+J)
   70       IF (JX.EQ.0) GO TO 120
*            U(1)=Q(JULEV+J)
*            U(2)=Q(JULEV+JX)
*            V(1)=Q(JVLEV+J)+(SIZD2*2.*FACHV)
*            V(2)=Q(JVLEV+JX)-(SIZD2*2.*FACHV)
            U(1)=Q(JULEV+J)
            U(2)=U(1)
            U(3)=Q(JULEV+JX)
            U(4)=Q(JULEV+JX)
            V(1)=Q(JVLEV+J)+(SIZD2*2.*FACHV)
            V(4)=Q(JVLEV+JX)-(SIZD2*2.*FACHV)
            V(2)=(V(1)+V(4))/2.
            V(3)=V(2)
            IF(J.EQ.1)HALF=V(4)-V(3)
            IF(NNPAR.EQ.3)THEN
               CALL IGPID(1,' ',J,' ')
            ENDIF
            CALL GDRAWV(V,U,4)
***
           IF(INTFLA.EQ.-1)THEN
            ARROWS=(SIZD2*2.*FACHV)/10.
            NPO=IQ(JTICK+J)
            FRA=(V(2)-V(1))/NPO
            DO 80  KJI=1,NPO
               VVV(KJI)=V(1)-(FRA/2.)+(FRA*KJI)
               UUU(KJI)=U(1)
   80       CONTINUE
            CALL ISFACI(2)
            CALL ISFAIS(1)
            DO 90  KJI=1,NPO
               SVVX(1)=VVV(KJI)-ARROWS
               SVVX(2)=VVV(KJI)-ARROWS
               SVVX(3)=VVV(KJI)+ARROWS
               SVVY(1)=UUU(KJI)+ARROWS
               SVVY(2)=UUU(KJI)-ARROWS
               SVVY(3)=UUU(KJI)
               IF(NNPAR.EQ.3)THEN
                  CALL IGPID(1,'Tree',IQ(JNAM1+J),' ')
                  CALL IGPID(2,'Arrow',KJI+1,' ')
               ENDIF
               CALL IFA(3,SVVX,SVVY)
   90       CONTINUE
***            CALL GDRAWV(U,V,2)
            NPO=IQ(JMYLL+J)
            FRA=(V(4)-V(3))/NPO
            DO 100 KJI=1,NPO
               VVV(KJI)=V(3)-(FRA/2.)+(FRA*KJI)
               UUU(KJI)=U(3)
  100       CONTINUE
            CALL ISFACI(3)
            CALL ISFAIS(1)
            IGREEN=IGREEN+1
            IORGO=0
            DO 110 KJI=NPO,1,-1
               IORGO=IORGO+1
               SVVX(1)=VVV(KJI)+ARROWS
               SVVX(2)=VVV(KJI)+ARROWS
               SVVX(3)=VVV(KJI)-ARROWS
               SVVY(1)=UUU(KJI)+ARROWS
               SVVY(2)=UUU(KJI)-ARROWS
               SVVY(3)=UUU(KJI)
               IF(NNPAR.EQ.3)THEN
                  CALL IGPID(1,'Tree',IQ(JNAM1+IGREEN),' ')
                  CALL IGPID(2,'Arrow',-IORGO,' ')
               ENDIF
               CALL IFA(3,SVVX,SVVY)
  110       CONTINUE
           ENDIF
***
            JX=IQ(JBRO1+JX)
            GO TO 70
  120       CONTINUE
         ELSE IF (IQ(JSCA1+J).EQ.-1) THEN
           IF(INTFLA.EQ.-1)THEN
            ARROWS=(SIZD2*2.*FACHV)/10.
            U(1)=Q(JULEV+J)
            U(2)=U(1)
            V(1)=Q(JVLEV+J)+(SIZD2*2.*FACHV)
            V(2)=V(1)+(SIZE*2.*FACHV)
            CALL GDRAWV(V,U,2)
            NPO=IQ(JTICK+J)
            FRA=(V(2)-V(1))/(NPO*2.)
            DO 130 KJI=1,NPO
               VVV(KJI)=V(1)-(FRA/2.)+(FRA*KJI)
               UUU(KJI)=U(1)
  130       CONTINUE
            CALL ISFACI(2)
            CALL ISFAIS(1)
            DO 140 KJI=1,NPO
               SVVX(1)=VVV(KJI)-ARROWS
               SVVX(2)=VVV(KJI)-ARROWS
               SVVX(3)=VVV(KJI)+ARROWS
               SVVY(1)=UUU(KJI)+ARROWS
               SVVY(2)=UUU(KJI)-ARROWS
               SVVY(3)=UUU(KJI)
               IF(NNPAR.EQ.3)THEN
                  CALL IGPID(1,'Tree',IQ(JNAM1+J),' ')
                  CALL IGPID(2,'Arrow',KJI+1,' ')
               ENDIF
               CALL IFA(3,SVVX,SVVY)
  140       CONTINUE
***            CALL GDRAWV(U,V,2)
           ENDIF
         ENDIF
  150 CONTINUE
C
C             Restore GDRAW calling parameters
C             and ZOOM internal parameters
C
  160 CONTINUE
      IF (ISEL4.NE.0) THEN
         GTHETA=SAVTHE
         GPHI=SAVPHI
         GU0=SAVU0
         GV0=SAVV0
         GSCU=SAVSCU
         GSCV=SAVSCV
         NGVIEW=0
         GZUA=SVGZUA
         GZVA=SVGZVA
         GZUB=SVGZUB
         GZVB=SVGZVB
         GZUC=SVGZUC
         GZVC=SVGZVC
      ENDIF
C
C             Restore original ranges previously saved
C
      PLTRNX=SAVPLX
      PLTRNY=SAVPLY
C
************      CALL IGRNG(PLTRNX,PLTRNY)
C
C             Delete tree structure on view bank 11
C
  170 CALL GDTR99(IVTREE)
      GO TO 999
C
  180 WRITE(CHMAIL,10000)KNAME
      CALL GMAIL(0,0)
10000 FORMAT(' GDTREE: VOLUME ',A4,' DOES NOT EXIST')
C
  999 END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR0(IVTREE,IROOT,LEVMAX,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Steering routine to build geometric structure tree       *
C.    *       (starting from JVOLUM structure until to have an         *
C.    *       optimized tree on IQ(JNAM1),IQ(JBRO1),etc.)              *
C.    *                                                                *
C.    *       IVTREE= view bank used to store tree structure (input)   *
C.    *       IROOT = root of the tree to be built up (input)          *
C.    *       IER   = error flag to detect nodes overflow (output)     *
C.    *                                                                *
C.    *    ==>Called by : GDFSPC, GDTREE                               *
C.    *       Author : P.Zanarini   *********                          *
C.    *       Revised: S.Giani      *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
      COMMON/GCLOCA/NLOCAL(2),JLOCA(19),JV
      COMMON /QUEST/ IQUEST(100)
      CHARACTER*4 NAMZ(16)
      SAVE NAMZ
      DATA NAMZ/'DNAM','DMOT','DXON','DBRO','DDUP','DSCA','DDVM','DPSM',
     +          'DNA1','DMO1','DXO1','DBR1','DDU1','DSC1','DUUU','DVVV'/
C.
C.    ------------------------------------------------------------------
C.
      IF (IVTREE.LE.NKVIEW.AND.LQ(JDRAW-IVTREE).NE.0) THEN
         Q(JDRAW+IVTREE)=1.
         CALL GDELET(IVTREE)
      ELSEIF (IVTREE.GT.NKVIEW) THEN
         NDL=IVTREE-NKVIEW
         CALL MZNEED(IXCONS,2*NDL,'G')
         IF(IQUEST(11).LT.0) GOTO 50
         CALL MZPUSH(IXCONS,JDRAW,NDL,NDL,'I')
         DO 10 I=NKVIEW+1,IVTREE-1
            Q(JDRAW+I)=1.
   10    CONTINUE
         NKVIEW=IVTREE
      ENDIF
C
C             Open view number IVTREE
C
      Q(JDRAW+IVTREE)=3.
      CALL MZBOOK(IXCONS,JV,JDRAW,-IVTREE,'DTRE',16,16,0,3,0)
      IF (JV.EQ.0) GO TO 40
*
      MAXNOD=16000
   20 ND=(10+MAXNOD)*16
      LLL=50+ND
      CALL MZNEED(IXCONS,LLL,'G')
      IQUEST(1)=0
      IF (IQUEST(11).LT.0) THEN
         MAXNOD=MAXNOD*0.9
         GO TO 20
      ENDIF
***   PRINT *,'NEW MAXNOD',MAXNOD
*
      DO 30 I=1,16
         CALL MZBOOK(IXCONS,LBANK,JV,-I,NAMZ(I),0,0,MAXNOD,3,0)
         JV=LQ(JDRAW-IVTREE)
         IF (LBANK.EQ.0) GO TO 40
   30 CONTINUE
C
      GO TO 60
C
   40 WRITE (CHMAIL,10000)
      CALL GMAIL(0,0)
      GO TO 999
C
   50 WRITE (CHMAIL,10100)
      CALL GMAIL(0,0)
      GO TO 999
C
   60 CONTINUE
C
C             1. Scan JVOLUM banks and build non-optimized tree structur
C                with arrays of pointers JNAM,JMOT,JXON,JBRO,JDUP,JDVM
C                Note : IQ(JDUP) is set >0 for multiplicity by positioni
C                       (normal case is 1); IQ(JDVM) is set <0 for
C                       multiplicity by dividing (or 0 if not)
C
C             2. Scan non-optimized tree and build optimized tree,
C                with arrays of pointers JNAM1,JMOT1,JXON1,JBRO1,JDUP1
C                Note : IQ(JDUP1) is set >0 for multiplicity by position
C                       or <0 for multiplicity by dividing (normal case
C
C
C             I = node pointer
C
C             IQ(JNAM+I) = node name (in non-optimized tree)
C             IQ(JMOT+I) = mother (or 0 if top of tree)
C             IQ(JXON+I) = first son (or 0 if any)
C             IQ(JBRO+I) = next brother (or 0 if any)
C             IQ(JDUP+I) = multiplicity (<0 for div. or >0 for pos.)
C             IQ(JSCA+I) = service word (scan flag, etc.)
C             IQ(JDVM+I) = division multiplicity (>0 for div. or 0 for n
C             IQ(JPSM+I) = position multiplicity (>1 for pos. or 1 for d
C
C             IQ(JNAM1+I) = node name (in optimized tree)
C             IQ(JMOT1+I) = mother (or 0 if top of tree)
C             IQ(JXON1+I) = first son (or 0 if any)
C             IQ(JBRO1+I) = next brother (or 0 if any)
C             IQ(JDUP1+I) = multiplicity (<0 for div. or >0 for pos.)
C             IQ(JSCA1+I) = service word (scan flag, etc.)
C             IQ(JULEV+I) = u-coordinate of the node in the tree
C             IQ(JVLEV+I) = v-coordinate of the node in the tree
C
      JNAM=LQ(JV-1)
      JMOT=LQ(JV-2)
      JXON=LQ(JV-3)
      JBRO=LQ(JV-4)
      JDUP=LQ(JV-5)
      JSCA=LQ(JV-6)
      JDVM=LQ(JV-7)
      JPSM=LQ(JV-8)
      JNAM1=LQ(JV-9)
      JMOT1=LQ(JV-10)
      JXON1=LQ(JV-11)
      JBRO1=LQ(JV-12)
      JDUP1=LQ(JV-13)
      JSCA1=LQ(JV-14)
      JULEV=LQ(JV-15)
      JVLEV=LQ(JV-16)
C
      DO 70  J=1,MAXNOD
         IQ(JSCA+J)=0
         IQ(JPSM+J)=1
         IQ(JDUP+J)=0
   70 CONTINUE
C
C             Scan JVOLUM banks and build non-optimized tree structure
C
C      WRITE (CHMAIL,5001)
C      CALL GMAIL(0,0)
      CALL GDTR2(IROOT,JMAX,LEVMAX,IER)
      IF (IER.NE.0) GO TO 999
      NUMNOD=JMAX
C
C             Scan non-optimized tree and build optimized tree
C
C      WRITE (CHMAIL,5002)
C      CALL GMAIL(0,0)
      CALL GDTR4
C      WRITE (CHMAIL,5003)
C      CALL GMAIL(0,0)
      CALL GDTR6(JMAX)
      NUMND1=JMAX
C
C      WRITE (CHMAIL,1500) NLEFT,MAXNOD,NUMNOD,NUMND1
C      CALL GMAIL(0,0)
C
10000 FORMAT (' ERROR IN GDTR0: MEMORY OVERFLOW IN CREATING A BANK')
10100 FORMAT (' ERROR IN GDTR0: MEMORY OVERFLOW IN PUSHING A BANK')
C
  999 END

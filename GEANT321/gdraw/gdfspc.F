*
* $Id$
*
* $Log$
* Revision 1.2  1998/09/25 09:23:48  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDFSPC(NAME,ISORT,INTER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a set of full specification of KNAME                *
C.    *       and all its sons (not duplicated);                       *
C.    *       if ISORT=1 plots will be alfabetically sorted;           *
C.    *       if INTER=1 it will prompt the user at each plot          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      LOGICAL MORE
      CHARACTER*4 CHKEY,NAME,NAME1,NAMOLD
C.
C.    ------------------------------------------------------------------
C.
C             Is NAME an existing volume ?
C
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF (IVO.LE.0) GO TO 999
C
C             Build tree structure using view bank 11
C
      IVTREE=NKVIEW+1
      CALL UCTOH (NAME, IROOT, 4, 4)
      CALL GDTR0 (IVTREE, IROOT, 0, IER)
      IF (IER.NE.0) GO TO 999
C
      JN=1
C
      DO 10 J=2,NUMND1
         IQ(JSCA1+J)=0
   10 CONTINUE
C
      IQ(JSCA1+1)=JN
      JMAX=1
      J=1
   20 CONTINUE
      JX=IQ(JXON1+IQ(JSCA1+J))
   30 IF (JX.EQ.0) GO TO 40
      JMAX=JMAX+1
      IQ(JSCA1+JMAX)=JX
      JX=IQ(JBRO1+JX)
      GO TO 30
   40 J=J+1
      IF (J.LE.JMAX) GO TO 20
C
      DO 50 I=1,JMAX
         IQ(JSCA1+I)=IQ(JNAM1+IQ(JSCA1+I))
   50 CONTINUE
C
      IF (ISORT.NE.1) GO TO 110
#if defined(CERNLIB_VAX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)||(defined(CERNLIB_LINUX) && !defined(CERNLIB_PPC))
C
C             Reverse bytes of word to be sorted
C
      CALL VXINVB(IQ(JSCA1+1),JMAX)
#endif
C
C             Perform bubble sort on names in ISCA1 (ascending order)
C
      NEXLAS=JMAX-1
      MORE=.TRUE.
C
   70 CONTINUE
      IF (.NOT.MORE) GO TO 90
      MORE=.FALSE.
      DO 80  I=1,NEXLAS
         IF (IQ(JSCA1+I).LE.IQ(JSCA1+I+1)) GO TO 80
         ITEMP=IQ(JSCA1+I)
         IQ(JSCA1+I)=IQ(JSCA1+I+1)
         IQ(JSCA1+I+1)=ITEMP
         MORE=.TRUE.
   80 CONTINUE
      GO TO 70
C
   90 CONTINUE
#if defined(CERNLIB_VAX)||defined(CERNLIB_MSDOS)||defined(CERNLIB_WINNT)||(defined(CERNLIB_LINUX) && !defined(CERNLIB_PPC))
C
C             Reverse bytes of words that have been sorted
C
      CALL VXINVB(IQ(JSCA1+1),JMAX)
#endif
C
C             Draw specs stored in ISCA1
C
  110 CONTINUE
      NAMOLD=' '
      DO 120 I=1,JMAX
         KNAME=IQ(JSCA1+I)
         CALL UHTOC(KNAME,4,NAME1,4)
         IF (NAME1.EQ.NAMOLD) GO TO 120
         NAMOLD=NAME1
         WRITE (CHMAIL,1000) NAME1
         CALL GMAIL(0,0)
         CHKEY=' '
         NCH=0
         IF (INTER.EQ.1) THEN
            CALL IGTERM
            CALL KUALFA
            CALL KUPROC('<CR>, NO, STOP',CHKEY,NCH)
            IF (CHKEY.EQ.'STOP') GO TO 999
         ENDIF
         CALL ICLRWK(0,0)
         IF (NCH.GT.0) GO TO 120
         CALL GDSPEC(NAME1)
  120 CONTINUE
C
C             Delete tree structure on view bank 11
C
      CALL GDTR99(IVTREE)
C
 1000 FORMAT (' DRAWING SPEC OF ',A4)
  999 RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GPRINT(CHNAME,NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print data structures                         *
C.    *                                                                *
C.    *       CHNAME   name of a data structure                        *
C.    *       NUMB     data structure number                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT  ,<GXINT> GINC4               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      CHARACTER*4 KNAMES(11),NAME
      DIMENSION NHAMES(11)
      CHARACTER*(*) CHNAME
      SAVE IFIRST,KNAMES,NHAMES
      DATA IFIRST/0/
      DATA KNAMES/'DIGI','JXYZ','HITS','KINE','MATE','VOLU'
     +           ,'ROTM','SETS','TMED','PART','VERT'/
C.
C.    ------------------------------------------------------------------
C.
      NAME=CHNAME
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH(KNAMES,NHAMES,4,44)
      ENDIF
C
      CALL GLOOK(NAME,NHAMES,11,JUMP)
C
      IF(JUMP.EQ. 1) CALL GPDIGI ('*','*')
      IF(JUMP.EQ. 2) CALL GPJXYZ (NUMB)
      IF(JUMP.EQ. 3) CALL GPHITS ('*','*')
      IF(JUMP.EQ. 4) CALL GPKINE (NUMB)
      IF(JUMP.EQ. 5) CALL GPMATE (NUMB)
      IF(JUMP.EQ. 6) CALL GPVOLU (NUMB)
      IF(JUMP.EQ. 7) CALL GPROTM (NUMB)
      IF(JUMP.EQ. 8) CALL GPSETS ('*','*')
      IF(JUMP.EQ. 9) CALL GPTMED (NUMB)
      IF(JUMP.EQ.10) CALL GPPART (NUMB)
      IF(JUMP.EQ.11) CALL GPVERT (NUMB)
C
      END

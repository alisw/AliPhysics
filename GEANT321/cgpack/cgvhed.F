*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGVHED(XYZ,IFACE,IPF,DFACE,ITREE,ISTACK,ALEFT,ARIGHT)
************************************************************************
*                                                                      *
*     Name: CGVHED                                                     *
*     Author: E. Chernyaev                       Date:    08.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Hidden line removal algorithm                          *
*               (tree search is used for min-max test)                 *
*                                                                      *
*     References: CGVEDG, CGVTEF                                       *
*                                                                      *
*     Input: XYZ(3,*) - nodes coordinates                              *
*         IFACE(*) - faces                                             *
*           IPF(*) - pointers to faces beginings                       *
*       DFACE(6,*) - faces min-max                                     *
*       ITREE(4,*) - tree of faces min-max                             *
*        ISTACK(*) - stack for tree search                             *
*         ALEFT(*) - min-max of left subtrees                          *
*        ARIGHT(*) - min-max of rigth subtrees                         *
*                                                                      *
*     Output: none                                                     *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgcedg.inc"
      REAL            XYZ(3,*),DFACE(6,*),ALEFT(*),ARIGHT(*)
*SG
      INTEGER         IPF(*)
      INTEGER         IFACE(*)
      INTEGER         ISTACK(*),ITREE(4,*)
      INTEGER         INDLFT(5),INDRGT(5)
*SG
      DATA            INDLFT/3,4,5,1,2/,INDRGT/3,4,1,2,1/
*-
      JS     = 0
      JT     = 1
*           T R E E   S E A R C H
  200 NF     = ITREE(1,JT)
      IND    = ITREE(4,JT)
      IF (DEDGE(IND) .GT. DFACE(IND,NF))      GOTO 300
      DO 210 I=1,5
        IF (DEDGE(I) .GT. DFACE(I,NF))        GOTO 220
  210   CONTINUE
      NIPF=IPF(NF)
      CALL CGVTEF(IFACE(NIPF),XYZ,IVIS)
      IF (IVIS .LT. 0)                        GOTO 900
*           S T O R E   L E F T   L I N K   I N   S T A C K  (IF NEED)
  220 IF (DEDGE(INDLFT(IND)) .GT. ALEFT(JT))  GOTO 300
      JS     = JS + 1
      ISTACK(JS) = ITREE(2,JT)
*           S T E P   T O   R I G H T
  300 IF (DEDGE(INDRGT(IND)) .GT. ARIGHT(JT)) GOTO 400
      JT     = ITREE(3,JT)
      GOTO 200
*           R E S T O R E   L E F T   L I N K   F R O M   S T A C K
  400 IF (JS .EQ. 0)                          GOTO 999
      JT     = ISTACK(JS)
      JS     = JS - 1
      GOTO 200
*
  900 NPART = -1
  999 RETURN
      END

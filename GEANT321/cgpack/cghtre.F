*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGHTRE(NFACE,DFACE,IORDER,ITREE,ALEFT,ARIGHT)
************************************************************************
*                                                                      *
*     Name: CGHDFA                                                     *
*     Author: E. Chernyaev                       Date:    07.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: build tree of faces min-max                            *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: NFACE      - number of faces                              *
*            DFACE(6,*) - min-max of faces                             *
*            IORDER     - work array                                   *
*                                                                      *
*     Output: ITREE(4,*) - tree of faces min-max                       *
*             ALEFT(*)   - min-max of left subtree                     *
*             ARIGHT(*)  - min-max of right subtree                    *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL            DFACE(6,*),ALEFT(*),ARIGHT(*)
      REAL            RNDM(1)
*SG
      INTEGER         IORDER(*),ITREE(4,*)
      INTEGER         INDLFT(5),INDRGT(5)
*SG
      DATA            INDLFT/3,4,5,1,2/,INDRGT/3,4,1,2,1/
*-
      DO 100 I=1,NFACE
        IORDER(I)  = I
  100   CONTINUE
*
**           T R E E   B U I L D
*
      K      = NFACE + 1
      IND    = 0
      JFREE  = 1
      DO 500 I=1,NFACE
        K      = K - 1
        CALL GRNDM(RNDM,1)
        IRNDM  = INT(RNDM(1)*K) + 1
        KF     = IORDER(IRNDM)
        IORDER(IRNDM) = IORDER(K)
        IF (I .EQ. 1)                   GOTO 400
        IT     = 1
  200   JT     = IT
        NF     = ITREE(1,JT)
        IND    = ITREE(4,JT)
        IF (DFACE(IND,KF) .GT. DFACE(IND,NF)) GOTO 300
*             S T E P   T O   L E F T
        INDL    = INDLFT(IND)
        IF (DFACE(INDL,KF) .GT. ALEFT(JT))    ALEFT(JT)=DFACE(INDL,KF)
        IT     = ITREE(2,JT)
        IF (IT .NE. 0)                  GOTO 200
        ITREE(2,JT) = JFREE
        GOTO 400
*             S T E P   T O   R I G H T
  300   INDR    = INDRGT(IND)
        IF (DFACE(INDR,KF) .GT. ARIGHT(JT))   ARIGHT(JT)=DFACE(INDR,KF)
        IT     = ITREE(3,JT)
        IF (IT .NE. 0)                  GOTO 200
        ITREE(3,JT) = JFREE
*             S E T   N E W   T R E E   N O D E
  400   IND    = IND + 1
        IF (IND .EQ. 6)                 IND = 1
        ITREE(1,JFREE) = KF
        ITREE(2,JFREE) = 0
        ITREE(3,JFREE) = 0
        ITREE(4,JFREE) = IND
        ALEFT (JFREE)  =-99999.
        ARIGHT(JFREE)  =-99999.
        JFREE  = JFREE+1
  500   CONTINUE
  999 RETURN
      END

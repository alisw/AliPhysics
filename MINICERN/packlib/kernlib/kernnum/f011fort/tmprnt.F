*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:04  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          TMPRNT(NAME,N,IDIM,K)
          CHARACTER*6         NAME
          LOGICAL             MFLAG,    RFLAG
          IF(NAME(2:2) .EQ. 'S') THEN
             CALL KERMTR('F012.1',LGFILE,MFLAG,RFLAG)
          ELSE
             CALL KERMTR('F011.1',LGFILE,MFLAG,RFLAG)
          ENDIF
          IF(NAME(3:6) .EQ. 'FEQN') ASSIGN 1002 TO IFMT
          IF(NAME(3:6) .NE. 'FEQN') ASSIGN 1001 TO IFMT
          IF(MFLAG) THEN
             IF(LGFILE .EQ. 0) THEN
                IF(NAME(3:6) .EQ. 'FEQN') THEN
                   WRITE(*,IFMT) NAME, N, IDIM, K
                ELSE
                   WRITE(*,IFMT) NAME, N, IDIM
                ENDIF
             ELSE
                IF(NAME(3:6) .EQ. 'FEQN') THEN
                   WRITE(LGFILE,IFMT) NAME, N, IDIM, K
                ELSE
                   WRITE(LGFILE,IFMT) NAME, N, IDIM
                ENDIF
             ENDIF
          ENDIF
          IF(.NOT. RFLAG) CALL ABEND
          RETURN
1001      FORMAT(7X, 31H PARAMETER ERROR IN SUBROUTINE , A6,
     +             27H ... (N.LT.1 OR IDIM.LT.N).,
     +             5X, 3HN =, I4, 5X, 6HIDIM =, I4, 1H. )
1002      FORMAT(7X, 31H PARAMETER ERROR IN SUBROUTINE , A6,
     +             37H ... (N.LT.1 OR IDIM.LT.N OR K.LT.1).,
     +             5X, 3HN =, I4, 5X, 6HIDIM =, I4, 5X, 3HK =, I4,1H.)
          END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:48:49  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
      SUBROUTINE F010PR(NAME,N,IDIM,K,KPRNT)
      CHARACTER*6 NAME
      LOGICAL MFLAG,RFLAG
C
C     ******************************************************************
C
C     PRINT ROUTINE FOR PARAMETER ERRORS IN MATRIX SUBROUTINES $EQINV,
C     $EQN, $INV (WHERE $ IS A LETTER SPECIFYING THE ARITHMETIC TYPE).
C
C     NAME         (CHARACTER*6) NAME OF THE CALLING ROUTINE.
C
C     N,IDIM,K     PARAMETERS OF THE CALLING ROUTINE (WITH K=0 IF K IS
C                  NOT TO BE PRINTED).
C
C     KPRNT        PRINT FLAG FOR K (K IS NOT PRINTED IF KPRNT=0).
C
C     ******************************************************************
C
C  START.
      CALL KERMTR('F010.1',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
         IF(LGFILE.EQ.0)  THEN
            IF(KPRNT.EQ.0) WRITE(*,2000) NAME,N,IDIM
            IF(KPRNT.NE.0) WRITE(*,2001) NAME,N,IDIM,K
         ELSE
            IF(KPRNT.EQ.0) WRITE(LGFILE,2000) NAME,N,IDIM
            IF(KPRNT.NE.0) WRITE(LGFILE,2001) NAME,N,IDIM,K
         ENDIF
      ENDIF
      IF(.NOT. RFLAG) CALL ABEND
      RETURN
C
 2000 FORMAT( 7X, 11HSUBROUTINE , A6, 14H ... PARAMETER,
     *        29H ERROR (N.LT.1 OR N.GT.IDIM).,
     *        6X, 3HN =, I4, 6X, 6HIDIM =, I4, 1H. )
 2001 FORMAT( 7X, 11HSUBROUTINE , A6, 14H ... PARAMETER,
     *        39H ERROR (N.LT.1 OR N.GT.IDIM OR K.LT.1).,
     *        6X, 3HN =, I4, 6X, 6HIDIM =, I4, 6X, 3HK =, I4, 1H. )
      END

*
* $Id$
*
* $Log$
* Revision 1.1  2005/07/29 16:15:46  jmb
* Include the various CERNLIB functions jimmy needs
*
* Revision 1.1.1.1  1996/04/01 15:02:52  mclareni
* Mathlib gen
*
      SUBROUTINE JMMLPT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT
      LOGICAL LMF,LRF
      SAVE

      IF(ERC(5:6).NE.'.0') THEN
        CALL JMMLMR(ERC,MLG,LMF,LRF)
      ELSE
        LMF=.TRUE.
        LRF=.FALSE.
      ENDIF
      IF(LMF) THEN
        LT=JMLNOC(TEXT)
        IF(MLG .LT. 1) WRITE(  *,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
        IF(MLG .GE. 1) WRITE(MLG,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
      ENDIF
      IF(.NOT.LRF) CALL JMABEND
      RETURN
100   FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END

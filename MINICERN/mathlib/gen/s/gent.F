*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/15 10:02:56  mclareni
* Change format statement to compile also on OSF, and add a date
*
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
      SUBROUTINE GENT ( LUN)
C.
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       PRINT PAM TITLE ON LUN                                   *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      WRITE (LUN,1000)
      RETURN
 1000 FORMAT('  GENT  COMPILED FROM PAM FILE= ',
     + ' GEN      1.48/00  150496' ,
     +/)
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GREND (LUN)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to close an RZ file                              *
C.    *                                                                *
C.    *       CHDIR    Top directory name of the file                  *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      CHARACTER*5 CHDIR
C.
C.    ------------------------------------------------------------------
C.
      IF(LUN.LT.10) THEN
         WRITE(CHDIR,'(''LUN'',I1)') LUN
      ELSE
         WRITE(CHDIR,'(''LUN'',I2)') LUN
      ENDIF
      CALL RZEND(CHDIR)
C
      END

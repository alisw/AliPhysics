*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GRMDIR(CHDIR,CHOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to create a "standard" GEANT subdirectory        *
C.    *                                                                *
C.    *       CHDIR Subdirectory name                                  *
C.    *                                                                *
C.    *       CHOPT is a character string which may be                 *
C.    *          ' ' To create a subdirectory                          *
C.    *          'S' To create a subdirectory and set the new          *
C.    *              Current Directory to this directory.              *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      CHARACTER*(*) CHDIR,CHOPT
      CHARACTER*8 CHTAGS(2)
      DATA CHTAGS/'Object','Version'/
C.
C.    ------------------------------------------------------------------
C.
      CALL UOPTC(CHOPT,'S',IOPTS)
*
      CALL RZMDIR(CHDIR,2,'HI',CHTAGS)
      IF(IOPTS.NE.0)THEN
         CALL RZCDIR(CHDIR,' ')
      ENDIF
*
      END

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
      SUBROUTINE GOPEN (LUN,LUNTYP,LEN,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to open I/O units                                *
C.    *                                                                *
C.    *       LUN      Logical unit number                             *
C.    *       LUNTYP   ZEBRA flag    I,O,X  (must be CHARACTER)        *
C.    *                'I' LUN is an INPUT file                        *
C.    *                'O' LUN is an OUTPUT file                       *
C.    *                'X' LUN is in exchange mode format              *
C.    *       LEN      Buffer length (see FZFILE)                      *
C.    *       IER      error flag                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
      CHARACTER*(*) LUNTYP
C.
C.    ------------------------------------------------------------------
C.
      IF (NUNITS.LT.5)THEN
         CALL FZFILE(LUN,LEN,LUNTYP)
         NUNITS = NUNITS + 1
         LUNITS(NUNITS) = LUN
         IF(INDEX(LUNTYP,'I').NE.0)LUNITS(NUNITS)=-LUN
         IER    = 0
      ELSE
         IER    = 1
      ENDIF
C
      END

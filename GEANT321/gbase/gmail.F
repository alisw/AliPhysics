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
      SUBROUTINE GMAIL(LINBEF,LINAFT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Send a mail on the output device (usually unit LOUT)     *
C.    *       The mail is in character array CHMAIL of /GCMAIL/        *
C.    *       Max length is 132 chars and trailing blanks are stripped *
C.    *       LINBEF lines are skipped before mail and LINAFT after    *
C.    *                                                                *
C.    *    ==>Called by : many routines                                *
C.    *       Authors : R.Brun, P.Zanarini   *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
C.
C.    ------------------------------------------------------------------
C.
      DO 10 I=1,LINBEF
        WRITE (LOUT,1000)
   10 CONTINUE
      DO 20 NCH=132,1,-1
        IF (CHMAIL(NCH:NCH).NE.' ') GO TO 30
   20 CONTINUE
      NCH=1
   30 CONTINUE
      WRITE (LOUT,2000) CHMAIL(1:NCH)
      DO 40 I=1,LINAFT
        WRITE (LOUT,1000)
   40 CONTINUE
 1000 FORMAT (1X)
 2000 FORMAT (A)
      END

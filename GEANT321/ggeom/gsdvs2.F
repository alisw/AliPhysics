*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSDVS2(CHNAME,CHMOTH,STEP,IAXIS,C0,NUMED)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    To be replaced in the future by GSDVT2                      *
C.    *                                                                *
C.    ******************************************************************
#include "geant321/gcunit.inc"
C.
      CHARACTER*4 CHNAME,CHMOTH
      WRITE(CHMAIL,
     +'('' *** GSDVS2: Obsolete Routine -- Call GSDVT2 instead'')')
      CALL GMAIL(1,0)
      WRITE(CHMAIL,
     +'(''            This code will be removed from next version'')')
      CALL GMAIL(0,1)
      CALL GSDVT2(CHNAME,CHMOTH,STEP,IAXIS,C0,NUMED,0)
      END

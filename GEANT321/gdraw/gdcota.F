*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_CG)
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCOTA
*
**********************************************************************
*                                                                    *
* RAY-TRACING                                                        *
* Function: Create colour's code table for ray-tracing and CG        *
*                                                                    *
* Author: S. Giani                                                   *
*                                                                    *
**********************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
*
      PARAMETER (NLSHAD=25,HINCRM=0.5/NLSHAD)
* red
      ID=16
      SL=0.2
      DO 10 I=1,NLSHAD
         ID=ID+1
         SL=SL+HINCRM
         CALL IGHTOR(0.,SL,.9,R,G,B)
         CALL ISCR(1,ID,R,G,B)
   10 CONTINUE
* yellow
      ID=41
      SL=0.2
      DO 20 I=1,NLSHAD
         ID=ID+1
         SL=SL+HINCRM
         CALL IGHTOR(60.,SL,.9,R,G,B)
         CALL ISCR(1,ID,R,G,B)
   20 CONTINUE
* green
      ID=66
      SL=0.2
      DO 30 I=1,NLSHAD
         ID=ID+1
         SL=SL+HINCRM
         CALL IGHTOR(120.,SL,.9,R,G,B)
         CALL ISCR(1,ID,R,G,B)
   30 CONTINUE
* light blue
      ID=91
      SL=0.2
      DO 40 I=1,NLSHAD
         ID=ID+1
         SL=SL+HINCRM
         CALL IGHTOR(180.,SL,.9,R,G,B)
         CALL ISCR(1,ID,R,G,B)
   40 CONTINUE
* blue
      ID=116
      SL=0.2
      DO 50 I=1,NLSHAD
         ID=ID+1
         SL=SL+HINCRM
         CALL IGHTOR(240.,SL,.9,R,G,B)
         CALL ISCR(1,ID,R,G,B)
   50 CONTINUE
* violet
      ID=141
      SL=0.2
      DO 60 I=1,NLSHAD
         ID=ID+1
         SL=SL+HINCRM
         CALL IGHTOR(300.,SL,.9,R,G,B)
         CALL ISCR(1,ID,R,G,B)
   60 CONTINUE
*
  999 END
#endif

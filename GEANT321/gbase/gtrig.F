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
      SUBROUTINE GTRIG
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Steering routine to process one event                    *
C.    *    ==>Called by : GRUN   ,<GXINT> GINC4                        *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvdma.inc"
#include "geant321/gcflag.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
C.
C.    ------------------------------------------------------------------
C.
C.            Kinematics
C.
      IF(RAYTRA.NE.1.)THEN
#if !defined(CERNLIB_USRJMP)
       CALL GUKINE
#endif
#if defined(CERNLIB_USRJMP)
       CALL JUMPT0(JUKINE)
#endif
      ELSE
       CALL GDKINE
      ENDIF
      IF (IEOTRI.NE.0)  GO TO  99
C.
C.            Tracking
C.
#if !defined(CERNLIB_USRJMP)
      CALL GUTREV
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUTREV)
#endif
      IF (RAYTRA.EQ.1.) GO TO  99
      IF (IEOTRI.NE.0)  GO TO  99
C.
C.            Digitisation
C.
#if !defined(CERNLIB_USRJMP)
      CALL GUDIGI
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUDIGI)
#endif
      IF (IEOTRI.NE.0)  GO TO  99
C.
C.            User end of event control routine
C.
#if !defined(CERNLIB_USRJMP)
      CALL GUOUT
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUOUT)
#endif
C
  99  RETURN
      END

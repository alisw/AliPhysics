*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GDKINE
*
************************************************************************
*                                                                      *
*            RAY-TRACING                                               *
*            Kinematics for the light rays.                            *
*            Parallel and perspective views.                           *
*                                                                      *
*            Author: S.Giani.                                          *
*                                                                      *
************************************************************************
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconst.inc"
********************************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcspee.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcrayt.inc"
********************************************************************************
      DIMENSION VERTEX(6),PLAB(3)
*
*     -----------------------------------------------------------------
*
      XMIP=0.
      CALL VZERO(VERTEX,6)
********************************************************************************
       CALL UCTOH('PERS',IPERS,4,4)
       XCOSXS=(SIN(GTHETA*DEGRAD))*(COS(GPHI*DEGRAD))
       YCOSYS=(SIN(GTHETA*DEGRAD))*(SIN(GPHI*DEGRAD))
       ZCOSZS=COS(GTHETA*DEGRAD)
*       print *,XCOSXS,YCOSYS,ZCOSZS,'direction from gukine'
       VDX=XCOSXS
       VDY=YCOSYS
       VDZ=ZCOSZS
       IK=48
       VERTEX(1)=XPINTS
       VERTEX(2)=YPINTS
       VERTEX(3)=ZPINTS
       IF(IPERS.EQ.IPRJ)THEN
        CONMOD=1./SQRT(((XPINTS-FPINTX)**2)+((YPINTS-FPINTY)**2)+
     +                 ((ZPINTS-FPINTZ)**2))
        XCOSXS=-(XPINTS-FPINTX)*CONMOD
        YCOSYS=-(YPINTS-FPINTY)*CONMOD
        ZCOSZS=-(ZPINTS-FPINTZ)*CONMOD
       ENDIF
       PLAB(1) = -PKINE(1)*XCOSXS
       PLAB(2) = -PKINE(1)*YCOSYS
       PLAB(3) = -PKINE(1)*ZCOSZS
********************************************************************************
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
      CALL GSKINE(PLAB,IK,NVERT,0,0,NT)
*
      END

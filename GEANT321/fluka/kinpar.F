*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE KINPAR.FOR
*COPY KINPAR
*
*=== kinpar ===========================================================*
*
      SUBROUTINE KINPAR(PXYZ,E,RCPMV)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      COMMON /FKKINP/AMTA,AMIT,PXTA,PYTA,PZTA,PSX,PSY,PSZ,PGES2,PGES,
     *PGESM,PTA2,PTAM,UMO2,CCX,CCY,CCZ,CXTA,CYTA,CZTA,GA,BGA
C*** CALCULATION OF KINEMATIC PARAMETERS
      DIMENSION PXYZ(10),E(2)
      PXTA=PXYZ(  2)
      PYTA=PXYZ(  6)
      PZTA=PXYZ( 10)
      PSX=PXYZ( 1)+P XTA
      PSY=PXYZ( 5)+P YTA
      PSZ=PXYZ( 9)+P ZTA
      PGES2=PSX*PSX+PSY*PSY+PSZ*PSZ
      PGES=SQRT(PGES2)
      PGESM=1.D0/PGES
      PTA2=PXTA*PXTA+PYTA*PYTA+PZTA*PZTA
      PTAM=1.D0/SQRT(PTA2)
      UMO2=((E(1)+E(2))**2-PGES2)/(RCPMV*RCPMV)
      UMO2=UMO2*0.000001D0
      CCX=PS X*PGESM
      CCY=PS Y*PGESM
      CCZ=PS Z*PGESM
      CXTA=P XTA*PTAM
      CYTA=P YTA*PTAM
      CZTA=P ZTA*PTAM
      GA=(E(1)+E(2))*0.001D0/(RCPMV*SQRT(UMO2))
      BGA=SQRT(GA*GA-1.D0)
      RETURN
      END

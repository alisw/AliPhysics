*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE DOST.FOR
*COPY DOST
*                                                                      *
*=== dost==============================================================*
*                                                                      *
      FUNCTION DOST(I,Z)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C--------------------------------------------------------------------
C SUBNAME = DOST
C--------------------------------------------------------------------
#include "geant321/eva0.inc"
C--------------------------------------------------------------------
      IF(Z-70.D0)1,2,2
    2 DOST=T(I,7)
    3 RETURN
    1 IF(Z-10.D0)5,5,6
    5 DOST=T(I,1)
      GOTO3
    6 N=.1D0*Z+1.D0
      X=10*N
      X=(X-Z)*.1D0
      DOST=X*T(I,N-1)+(1.D0-X)*T(I,N)
      GOTO3
      END

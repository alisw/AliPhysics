*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE ENRG.FOR
*COPY ENRG
*                                                                      *
*=== enrg =============================================================*
*                                                                      *
      FUNCTION ENRG(A,Z)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C--------------------------------------------------------------------
C SUBNAME = ENRG
C--------------------------------------------------------------------
      PARAMETER ( O16OLD = 931.145  D+00 )
      PARAMETER ( O16NEW = 931.19826D+00 )
      PARAMETER ( O16RAT = O16NEW / O16OLD )
      PARAMETER ( C12NEW = 931.49432D+00 )
      PARAMETER ( ADJUST = -8.322737768178909D-02 )
#include "geant321/eva0.inc"
      LOGICAL LFIRST
      SAVE LFIRST, EXHYDR, EXNEUT
      DATA LFIRST / .TRUE. /
      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         EXHYDR = FKENER ( ONEONE, ONEONE )
         EXNEUT = FKENER ( ONEONE, ZERZER )
      END IF
C--------------------------------------------------------------------
  100 IZ0 = NINT (Z)
      IF ( IZ0 .LE. 0 ) THEN
         ENRG = A * EXNEUT
         RETURN
      END IF
      N   = NINT (A-Z)
      IF ( N .LE. 0 ) THEN
         ENRG = Z * EXHYDR
         RETURN
      END IF
      AM2ZOA= (A-Z-Z)/A
      AM2ZOA=AM2ZOA*AM2ZOA
      A13 = RMASS(NINT(A))
      AM13 = 1.D0/A13
      EV=-17.0354D0*(1.D0 -1.84619D0*AM2ZOA)*A
      ES= 25.8357D0*(1.D0 -1.712185D0*AM2ZOA)*
     &    (1.D0 -0.62025D0*AM13*AM13)*
     &    (A13*A13 -.62025D0)
      EC= 0.799D0*Z*(Z-1.D0)*AM13*(((1.5772D0*AM13 +1.2273D0)*
     &    AM13-1.5849D0)*
     &    AM13*AM13 +1.D0)
      EEX= -0.4323D0*AM13*Z**1.3333333D0*
     &   (((0.49597D0*AM13 -0.14518D0)*AM13 -0.57811D0) * AM13 + 1.D0)
      ENRG  =8.367D0*A -0.783D0*Z +EV +ES +EC +EEX +CAM2(IZ0)+CAM3(N)
      ENRG  = ( ENRG + A * O16OLD ) * O16RAT - A * ( C12NEW - ADJUST )
      ENRG  = MIN ( ENRG, Z * EXHYDR + ( A - Z ) * EXNEUT )
      RETURN
      END

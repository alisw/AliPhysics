*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE XINNEU.FOR
*COPY XINNEU
*
*=== xinneu ===========================================================*
*
      FUNCTION XINNEU ( E, ZTAR, BBETAA )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/isotop.inc"
#include "geant321/xsepar.inc"
*
      IZ = NINT ( ZTAR )
      EEE = 0.015D+00 * ( 1.D+00 + 1.D-04 * ZTAR * ZTAR )
      TMPZZN = 2.D-01 * ZZNXSE (IZ)
      EEE = MAX ( EEE, TMPZZN )
      IF ( E .GT. 0.1D+00 ) THEN
         ATAR = 0.D+00
         DO 25 IS = ISONDX (1,IZ), ISONDX (2,IZ)
             ATAR = ATAR +  ISOMNM (IS) * ABUISO (IS)
   25    CONTINUE
         EMEV = E * 1.D+03
         XINNEU = SITSAO ( EMEV, IZ, ATAR )
      ELSE IF ( E .LT. EEE ) THEN
         X = 1.D+01 * EEE
         XSENE1 = AANXSE (IZ) * ( X - ZZNXSE (IZ) ) / X * ( 1.D+00
     &          + BBNXSE (IZ) * X / ( 1.D+00 + (CCNXSE(IZ)*X)**2 )
     &          + DDNXSE (IZ) * X * X / ( 1.D+00 + (EENXSE(IZ)*X)**4 ))
         BBCOFF = XSENE1 / ( 1.D+00 + BBETAA / EEE )
         AACOFF = BBCOFF * BBETAA
         XINNEU = AACOFF / E + BBCOFF
      ELSE IF ( E .GT. 0.012D+00 ) THEN
         X  = 1.D+01 * E
         XINNEU = AANXSE (IZ) * ( X - ZZNXSE (IZ) ) / X * ( 1.D+00
     &          + BBNXSE (IZ) * X / ( 1.D+00 + (CCNXSE(IZ)*X)**2 )
     &          + DDNXSE (IZ) * X * X / ( 1.D+00 + (EENXSE(IZ)*X)**4 ))
      ELSE
         XINNEU = 0.D+00
      END IF
      XINNEU = MAX ( XINNEU, ZERZER )
*=== End of function xinneu ===========================================*
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.42  by  S.Giani
*-- Author :
*$ CREATE FPOWER.FOR
*COPY FPOWER
*
*=== fpower ===========================================================*
*
      FUNCTION FPOWER ( IPOWER, NU, ANUAV )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on   05  may  1992   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 10-jun-92     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( TWOLOG = 0.6931471805599453D+00 )
*
      IF ( IPOWER .NE. 13 ) STOP 'FPOWER'
         IF ( ANUAV .LT. 2.D+00 ) THEN
            AKAPP0 = -0.93699278D+00*ANUAV**2 + 3.9082338D+00
         ELSE
            AKAPP0 = 4.2019D+00 / ANUAV**2.2835973D+00
         END IF
         IF ( NU .LE. 1 ) THEN
            FPOWER = AKAPP0
         ELSE
            AFACT  = ( NU - 1.D+00 ) / ( ANUAV + NU - 1.D+00 )
*  |  This correction (1 card) has been put by A.F. on 29-12-92, trying
*  |  to widen a bit the residual nuclei distribution without
*  |  spoiling the good agreement with the grey particle distribution
            AFACT  = AFACT / SQRT (ANUAV)
            AKAPPA = ( 1.D+00 - AFACT ) * AKAPP0 + AFACT
            APOWER = ANUAV / NU + 1.D+00 * ANUAV / ( ANUAV + 1.D+00 )
            DDNU   = NU
            FPOWER = APOWER + LOG ( AKAPPA ) / LOG ( DDNU )
         END IF
*=== End of function fpower ===========================================*
      RETURN
      END

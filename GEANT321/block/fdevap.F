*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:35  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE BDEVAP.FOR
*COPY BDEVAP
*
*=== bdevap ===========================================================*
*
      SUBROUTINE FDEVAP
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Block Data for the EVAPoration routines:                         *
*                                                                      *
*     Created on    20 may 1990    by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Modified from the original version of J.M.Zazula                 *
*     and, for cookcm, from a LAHET block data kindly provided by      *
*     R.E.Prael-LANL                                                   *
*                                                                      *
*     Last change on  28-jan-93    by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/eva0.inc"
#include "geant321/hettp.inc"
#include "geant321/hetc7.inc"
#include "geant321/inpflg.inc"
*
      DATA B01   / 8.D+00 /, Y01 / 1.5D+00 /
      DATA IANG1 / 1 /, IFISS1 / 0 /,  IB01 / 2 /, IGEOM1 / 0 /
      DATA ISTRA1 /0/, KEYDK1 /0/
      DATA NBERT1 /LUNBER/
      DATA COSTH1 /1.D+00/, SINTH1 /0.D+00/, COSPH1 /1.D+00/,
     &     SINPH1/0.D+00/
*
      B0=B01
      Y0=Y01
      IANG=IANG1
      IFISS=IFISS1
      IB0=IB01
      IGEOM=IGEOM1
      ISTRAG=ISTRA1
      KEYDK=KEYDK1
      NBERTP=NBERT1
      COSTH=COSTH1
      SINTH=SINTH1
      COSPHI=COSPH1
      SINPHI=SINPH1
      END

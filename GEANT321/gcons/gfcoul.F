*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      FUNCTION GFCOUL(Z)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Compute Coulomb correction for pair production and Brem       *
C.    *  REFERENCE : EGS MANUAL SLAC 210 - UC32 - JUNE 78              *
C.    *                        FORMULA 2.7.17                          *
C.    *                                                                *
C.    *    ==>Called by : GSMIXT                                       *
C.    *       Author    M.Maire *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      PARAMETER (ALPHA = 7.29927E-03)
C
      AZ2   = (ALPHA*Z)**2
      AZ4   =   AZ2 * AZ2
      FP    = ( 0.0083*AZ4 + 0.20206 + 1./(1.+AZ2) ) * AZ2
      FM    = ( 0.0020*AZ4 + 0.0369  ) * AZ4
      GFCOUL = FP - FM
      END

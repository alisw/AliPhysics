*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUINME(X,SH,P,IYES)
C.
C.    **********************************************
C.    *                                            *
C.    *    USER ROUTINE TO PROVIDE GINME FUNCTION  *
C.    *    FOR ALL USER SHAPES IDENTIFIED BY THE   *
C.    *    SHAPE NUMBER SH. POINT IS GIVEN IN X    *
C.    *    THE PARAMETERS ARE GIVEN IN P. IYES IS  *
C.    *    RETURNED 1 IF POINT IS IN, 0 IF POINT   *
C.    *    IS OUT AND LESS THAN ZERO IF SHAPE      *
C.    *    NUMBER IS NOT SUPPORTED.                *
C.    *                                            *
C.    *    ==>Called by : GINME                    *
C.    *                                            *
C.    **********************************************
C.
      DIMENSION X(3),P(*)
C
      IYES=-1
      END

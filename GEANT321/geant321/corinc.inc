*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CORINC_INC
#define CERNLIB_GEANT321_CORINC_INC
*
*
* corinc.inc
*
*                                                                      *
*=== corinc ===========================================================*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Modified by A. Ferrari: be sure that the parameter Inmax is >=   *
*     of the same parameter in the function Nudist                     *
*----------------------------------------------------------------------*
*
      PARAMETER ( INMAX = 20 )
      COMMON /FKCORI/ XSEA (INMAX),  XASEA (INMAX), FRAINC, ANUAV,
     &                IJTARG (INMAX+1), NSEA
 

#endif

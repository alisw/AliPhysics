*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_BAMJCM_INC
#define CERNLIB_GEANT321_BAMJCM_INC
*
*
* bamjcm.inc
*
*
*=== bamjcm ===========================================================*
*
 
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Bamjcm by A. Ferrari                                             *
*                           included in:                               *
*                                        bamjet                        *
*     When changing kmxjcm dimension look also at verein!!!!!          *
*----------------------------------------------------------------------*
*                                                                      *
      PARAMETER (KMXJCM = 100)
      COMMON /FKBAMJ/  RPX(KMXJCM),  RPY(KMXJCM),   RE(KMXJCM),
     &                KFR1(KMXJCM), KFR2(KMXJCM),   IV(KMXJCM)
 
 

#endif

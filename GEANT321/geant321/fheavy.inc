*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:38  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_FHEAVY_INC
#define CERNLIB_GEANT321_FHEAVY_INC
*
*
* fheavy.inc
*
*
*=== fheavy ===========================================================*
*
*----------------------------------------------------------------------*
*     Last change on 9 March 1993     by Alfredo Ferrari, INFN Milan   *
*     include file: fheavy copy       created 5/4/90 by A. Ferrari     *
*     changes: none                                                    *
*     included in the following subroutines or functions: not updated  *
*                                                                      *
*     description of the common block(s) and variable(s)               *
*                                                                      *
*     /fheavy/ is the storage for heavy secondaries created in the     *
*              nuclear evaporation                                     *
*        npheav     = number of secondaries                            *
*        kheavy(ip) = type of the secondary ip                         *
*                   ( 3 = deuteron, 4 = 3-H, 5 = 3-He, 6 = 4-He,       *
*                     7-12 = "Heavy" fragment specified by Ibheav and  *
*                     Icheav )                                         *
*        cxheav(ip) = direction cosine of the secondary ip             *
*                     with respect to x-axis                           *
*        cyheav(ip) = direction cosine of the secondary ip             *
*                     with respect to y-axis                           *
*        czheav(ip) = direction cosine of the secondary ip             *
*                     with respect to z-axis                           *
*        tkheav(ip) = kinetic energy of secondary ip                   *
*        pheavy(ip) = momentum of the secondary ip                     *
*        wheavy(ip) = weight of the secondary ip                       *
*        amheav(ip) = atomic masses of the six types of evaporated     *
*                     particles                                        *
*        anheav(kp) = name of the kp-type heavy particle               *
*        icheav(kp) = charge of the kp-type heavy particle             *
*        ibheav(kp) = mass number of the kp-type heavy particle        *
*   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   *
*   !!! There is now the possibility to produce up to 6 "heavy" !!!!   *
*   !!! fragments besides the residual nucleus recorded in      !!!!   *
*   !!! Resnuc: they are identified by indeces 7-12, of course  !!!!   *
*   !!! the corresponding physical properties (Z,A,m..) must be !!!!   *
*   !!! updated every time they are produced                    !!!!   *
*   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   *
*----------------------------------------------------------------------*
*
      PARAMETER ( MXHEAV = 100 )
      CHARACTER*8 ANHEAV
      COMMON / FKHEAV / CXHEAV (MXHEAV), CYHEAV (MXHEAV),
     &                  CZHEAV (MXHEAV), TKHEAV (MXHEAV),
     &                  PHEAVY (MXHEAV), WHEAVY (MXHEAV),
     &                  AMHEAV  ( 12 ) ,
     &                  KHEAVY (MXHEAV), ICHEAV  ( 12 ) ,
     &                  IBHEAV  ( 12 ) , NPHEAV
      COMMON / FKHEAC / ANHEAV  ( 12 )
 

#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_REAC_INC
#define CERNLIB_GEANT321_REAC_INC
*
*
* reac.inc
*
*
*=== reac =============================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 10 december 1991  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 10-dec-91     by    Alfredo Ferrari               *
*                                                                      *
*     This is the original common reac of Hadrin                       *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                CALUMV                                                *
*                DATESH                                                *
*                HADRIV                                                *
*                RCHANV                                                *
*                BLKDT3                                                *
*                HADRIN                                                *
*                IEFUN                                                 *
*                SIGINT                                                *
*                CALUMO                                                *
*                RCHANW                                                *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON /FKREAC/ UMO   ( 296), PLABF ( 296), SIIN  ( 296),
     +              WK    (5184), NRK  (2,268), NURE  (30,2)
 

#endif

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_RESNUC_INC
#define CERNLIB_GEANT321_RESNUC_INC
*
*
* resnuc.inc
*
*                                                                     *
*=== resnuc ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*     Include file Resnuc                                             *
*                                                                     *
*     Created on 20 april 1990  by            Alfredo Ferrari         *
*                                               INFN Milan            *
*                                                                     *
*     Last change on  9  march  1993                                  *
*                                                                     *
*     Included in the following routines:                             *
*                                                                     *
*            BDNOPT                                                   *
*            COREVT                                                   *
*            CORRIN                                                   *
*            DISTNU                                                   *
*            DRES                                                     *
*            EVDEEX                                                   *
*            EVENTQ                                                   *
*            EVENTV                                                   *
*            EVEVAP                                                   *
*            FEREVV                                                   *
*            FERHAV                                                   *
*            LOWRES                                                   *
*            NUCEVT                                                   *
*            NUCEVV                                                   *
*            NUCRIN                                                   *
*            NUCRIV                                                   *
*                                                                     *
*---------------------------------------------------------------------*
*
      LOGICAL LRNFSS
      COMMON /FKRESN/  AMNTAR, AMMTAR, AMNZM1, AMMZM1, AMNNM1, AMMNM1,
     &                   ANOW,   ZNOW, ANCOLL, ZNCOLL,   ERES,  EKRES,
     &                 AMNRES, AMMRES,  PTRES,  PXRES,  PYRES,  PZRES,
     &                 PTRES2,  KTARP,  KTARN, IGREYP, IGREYN,  ICRES,
     &                  IBRES, IEVAPL, IEVAPH, IEVNEU, IEVPRO, IEVDEU,
     &                 IEVTRI, IEV3HE, IEV4HE, IDEEXG,  IBTAR, ICHTAR,
     &                 IOTHER, LRNFSS
 

#endif

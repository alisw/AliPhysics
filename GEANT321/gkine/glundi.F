*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GLUNDI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        GEANT INTERFACE TO THE LUND MONTE-CARLO                 *
C.    *            (Initialization routine)                            *
C.    *                                                                *
C.    *  FOR TECHNICAL DETAILS OF THE LUND MONTE CARLO,SEE             *
C.    *  1. THE LUND MONTE CARLO FOR JET FRAGMENTATION - T.SJOSTRAND   *
C.    *     LU TP 82-3  and COMPUTER PHYSICS COMMUNICATION 27(1982)243 *
C.    *  2. THE LUND MONTE CARLO FOR E+E- JET PHYSICS  - T.SJOSTRAND   *
C.    *     LU TP 82-7 and COMPUTER PHYSICS COMMUNICATION 28(1983)229  *
C.    *                                                                *
C.    *  THE PARTICLES FORWARDED TO GEANT ARE                          *
C.    *  GAMMA,ELECTRON,NEUTRINO(E,MU,TAU),MU-,PI+,K+,K0S,K0L,         *
C.    *  PROTON,NEUTRON,SIG(+,-),XI(0,-),LAM,OM-                       *
C.    *  ....  AND THEIR ANTIPARTICLES                                 *
C.    *                                                                *
C.    *                                                                *
C.    *  THE USER MUST SUPPLY TWO PARAMETERS FOR EVENT GENERATION      *
C.    *  VIA     CALL LUEEVT(IFLUND,ECLUND)     WHERE                  *
C.    *  IFLUND=THE FLAVOUR CODE (1=U,2=D,3=S,4=C,5=B,6=T,0=MIX)       *
C.    *  ECLUND=THE TOTAL CMS ENERGY                                   *
C.    *                                                                *
C.    *     GLUNDI is automatically called by GLUND.                   *
C.    *     GLUNDI can be called by the user in the initialization     *
C.    *     routine UGINIT just before reading data cards (GFFGO)      *
C.    *     In this case the default parameters can be changed by      *
C.    *     the data cards KTYP,PMAS,PWID,MSTE,LUND and IDB            *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUKINE                               *
C.    *       Author    R.Hemingway , D.Ward *********                 *
C.    *       modified by R.Brun 18/11/85 ******                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcflag.inc"
#include "geant321/ludat1.inc"
#include "geant321/ludat3.inc"
#include "geant321/gclund.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF (IFINIT(3).EQ.0)THEN
*
*              Set defaults
*
         IFINIT(3) = 1
         IFLUND = 0
         ECLUND = 92.25
*
*--  Initial state radiative effects     MSTJ(107) =  1
*
         MSTJ(107)= 1
*
*--  Set some particles to be stable in LUND so that GEANT can handle
*--  their tracking and decays. (K0s, sigma+-, xi0-, lambda0, omega-)
*
         MDCY(PYCOMP(310) ,1)=0
         MDCY(PYCOMP(3222),1)=0
         MDCY(PYCOMP(3112),1)=0
         MDCY(PYCOMP(3322),1)=0
         MDCY(PYCOMP(3312),1)=0
         MDCY(PYCOMP(3122),1)=0
         MDCY(PYCOMP(3334),1)=0
*
         CALL FFKEY('LUND',IFLUND,2,'MIXED')
      ENDIF
C
      END

*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 09/01/95  16.20.45  by  S.Ravndal
*-- Author :
      SUBROUTINE GLUDKY
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        GEANT INTERFACE TO THE LUND MONTE-CARLO                 *
C.    *        FOR TAU DECAY GENERATION ONLY                           *
C.    *                                                                *
C.    *  THIS INTERFACE TO LUND                                        *
C.    *  A) INHIBITS DECAYS OF K0S,SIG+,SIG-,XI0,XI-,LAM,OM-           *
C.    *  B) ACCEPTS STABLE PARTICLES ONLY                              *
C.    *  C) CONVERTS THE PARTICLE CODES TO GEANT CONVENTIONS           *
C.    *                                                                *
C.    *  THUS THE PARTICLES FORWARDED TO GEANT ARE                     *
C.    *  GAMMA,ELECTRON,NEUTRINO(E,MU,TAU),MU-,PI+,K+,K0S,K0L,         *
C.    *  PROTON,NEUTRON,SIG(+,-),XI(0,-),LAM,OM-                       *
C.    *  ....  AND THEIR ANTIPARTICLES                                 *
C.    *                                                                *
C.    *       That routine can be called from the user routine GUDCAY  *
C.    *                                                                *
C.    *       Author    R.Hemingway  26 SEPT 1985                      *
C.    *                                                                *
C.    ******************************************************************
C.
C
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconsp.inc"
#include "geant321/lujets.inc"
#include "geant321/ludat3.inc"
      PARAMETER (NSEL=44)
      DIMENSION IPLUND(NSEL)
      DATA IPLUND/
     +    22,   -11,    11,    12,   -13,    13,   111,   211,
     +  -211,   130,   321,  -321,  2112,  2212, -2212,   310,
     +   221,  3122,  3222,  3212,  3112,  3322,  3312,  3334,
     + -2112, -3122, -3112, -3212, -3222, -3322, -3312, -3334,
     +   -15,    15,   411,  -411,   421,  -421,   431,  -431,
     +  4122,    24,   -24,    23/
*.
*.    ------------------------------------------------------------------
*.
      IF (IFINIT(6).EQ.0)THEN
*
*--  Set some particles to be stable in LUND so that GEANT can handle
*--  their tracking and decays. (K0s, sigma+-, xi0-, lambda0, omega-)
*
         MDCY(LUCOMP(310) ,1)=0
         MDCY(LUCOMP(3222),1)=0
         MDCY(LUCOMP(3112),1)=0
         MDCY(LUCOMP(3322),1)=0
         MDCY(LUCOMP(3312),1)=0
         MDCY(LUCOMP(3122),1)=0
         MDCY(LUCOMP(3334),1)=0
         IFINIT(6) = 1
      ENDIF
C
      IF(IPART.EQ.33.OR.IPART.EQ.34) THEN
         N=0
         EDKY=GETOT
         THDKY=ACOS(VECT(6))
         PHDKY=ATAN2(VECT(5),VECT(4))
*
*            for TAU+ TAU-
*            could be extended to all charm particles
*
         CALL LU1ENT(1,IPLUND(IPART),EDKY,THDKY,PHDKY)
         CALL LUEXEC
         NGKINE=0
         DO 20 I = 1,N
            IF(K(I,1).GE.1.AND.K(I,1).LE.10) THEN
               DO 10 IGE=1, NSEL
                  IF(IPLUND(IGE).EQ.K(I,2)) THEN
                     NGKINE=NGKINE+1
                     GKIN(1,NGKINE)=P(I,1)
                     GKIN(2,NGKINE)=P(I,2)
                     GKIN(3,NGKINE)=P(I,3)
                     GKIN(4,NGKINE)=P(I,4)
                     GKIN(5,NGKINE)=IGE
                     GPOS(1,NGKINE)=V(I,1)
                     GPOS(2,NGKINE)=V(I,2)
                     GPOS(3,NGKINE)=V(I,3)
                     TOFD(NGKINE)  =V(I,4)/(10.*CLIGHT)
                  ENDIF
   10          CONTINUE
            ENDIF
   20    CONTINUE
      ENDIF
*
      END

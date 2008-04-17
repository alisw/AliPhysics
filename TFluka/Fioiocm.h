#ifndef FIOIOCM_H
#define FIOIOCM_H

#include "Rtypes.h"
#include "cfortran.h"
extern "C" {

/*
*$ CREATE IOIOCM.ADD
*COPY IOIOCM
*
*=== Ioiocm ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2002-2006      by  Alfredo Ferrari                 *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     IOn-IOn CoMmon:                                                  *
*                                                                      *
*     Ion-Ion collision common for Fluka9x/Fluka200x....:              *
*                                                                      *
*     Last change  on  19-apr-06   by  Alfredo Ferrari, INFN-Milan     *
*                                                                      *
*     Description of the variable(s):                                  *
*                                                                      *
*        Eknion = laboratory kinetic  energy per nucleon  (GeV/amu)    *
*        Etnion = laboratory total    energy per nucleon  (GeV/amu)    *
*        Plnion = laboratory momentum        per nucleon  (GeV/c/amu)  *
*        Eexion = excitation energy of the projectile ion (GeV)        *
*        T12ion = half life of the projectile ion (for Eexion=0) (s)   *
*     Matprj(i) = list   of materials used as projectiles              *
*        Nmatpr = number of materials defined inside Matprj            *
*        Iproa  = the projectile mass   number                         *
*        Iproz  = the projectile proton number                         *
*        Iprom  = the projectile isomer number                         *
*                                                                      *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*     !!!! Note that the units are GeV/amu --> per unit mass  !!!!     *
*     !!!! with mass measured in amu (1 amu = Amuc12 GeV)     !!!!     *
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON / IOIOCM / EKNION, ETNION, PLNION, EEXION, T12ION,
     &                  MATPRJ (MXXMDF), NMATPR, IPROA , IPROZ , IPROM
      SAVE / IOIOCM /
*/
    typedef struct {
	Double_t eknion;
	Double_t etnion;
	Double_t plnion;
	Double_t eexion;
	Double_t t12ion;
	Int_t    matprj[mxxmdf];
	Int_t    nmatpr;
	Int_t    iproa;
	Int_t    iproz;
	Int_t    iprom;
    } ioiocmCommon;
#define IOIOCM COMMON_BLOCK(IOIOCM,ioiocm)
COMMON_BLOCK_DEF(ioiocmCommon, IOIOCM);

}
#endif

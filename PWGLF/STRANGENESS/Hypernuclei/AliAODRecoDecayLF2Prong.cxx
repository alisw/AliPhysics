/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed2-prong decay
// strongly based on AliAODRecoDecayHF2Prong  
// Author: Ramona Lea (ramona.lea@cern.ch)
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include "AliAODRecoDecayLF.h"
#include "AliAODRecoDecayLF2Prong.h"

ClassImp(AliAODRecoDecayLF2Prong)

//--------------------------------------------------------------------------
AliAODRecoDecayLF2Prong::AliAODRecoDecayLF2Prong() :
  AliAODRecoDecayLF()
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayLF2Prong::AliAODRecoDecayLF2Prong(AliAODVertex *vtx2,
						 Double_t *px,Double_t *py,Double_t *pz,
						 Double_t *d0,Double_t *d0err,Float_t dca) :
  AliAODRecoDecayLF(vtx2,2,0,px,py,pz,d0,d0err)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //
  SetDCA(dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayLF2Prong::AliAODRecoDecayLF2Prong(AliAODVertex *vtx2,
						 Double_t *d0,Double_t *d0err,Float_t dca) :
  AliAODRecoDecayLF(vtx2,2,0,d0,d0err)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //
  SetDCA(dca);
}
//--------------------------------------------------------------------------
AliAODRecoDecayLF2Prong::AliAODRecoDecayLF2Prong(const AliAODRecoDecayLF2Prong &source) :
  AliAODRecoDecayLF(source)
{
  //
  // Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecayLF2Prong &AliAODRecoDecayLF2Prong::operator=(const AliAODRecoDecayLF2Prong &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayLF::operator=(source);

  return *this;
}
